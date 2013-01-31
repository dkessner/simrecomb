//
// Population.cpp
//
// Copyright 2012 Darren Kessner
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
//

#include "Population.hpp"
#include "Random.hpp"
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iterator>
#include <set>
#include <fstream>


using namespace std;
using boost::shared_ptr;


//
// MatingDistribution
//


namespace {
struct HasLowerWeight
{
    bool operator()(const MatingDistribution::Entry& a, const MatingDistribution::Entry& b) 
    {
        return a.cumulativeWeight < b.cumulativeWeight;
    }
};
} // namespace


void MatingDistribution::push_back(double weight, const IndexPair& indexPair)
{
    totalWeight_ += weight;
    entries_.push_back(Entry(totalWeight_, indexPair));
}


const MatingDistribution::IndexPair& MatingDistribution::random_index_pair(const Random& random) const
{
    double roll = random.uniform(0, totalWeight_);

    vector<Entry>::const_iterator it = lower_bound(entries_.begin(), entries_.end(),
                                                   Entry(roll, make_pair(0,0)), HasLowerWeight());
    if (it == entries_.end())
    {
        cout << "totalWeight_: " << totalWeight_ << endl;
        cout << "roll: " << roll << endl;
        throw runtime_error("[MatingDistribution::Impl::random_index_pair()] This isn't happening.");
    }
    
    return it->indexPair;
}


bool operator==(const MatingDistribution::Entry& a, const MatingDistribution::Entry& b)
{
    return a.cumulativeWeight == b.cumulativeWeight &&
           a.indexPair.first == b.indexPair.first &&
           a.indexPair.second == b.indexPair.second;
}


bool operator!=(const MatingDistribution::Entry& a, const MatingDistribution::Entry& b)
{
    return !(a==b);
}


bool operator==(const MatingDistribution& a, const MatingDistribution& b)
{
    if (a.entries().size() != b.entries().size()) return false;

    for (size_t i=0; i<a.entries().size(); ++i)
        if (a.entries()[i] != b.entries()[i])
            return false;

    return true;
}


bool operator!=(const MatingDistribution& a, const MatingDistribution& b)
{
    return !(a==b);
}


ostream& operator<<(ostream& os, const MatingDistribution& md)
{
    // {<.6|0,0><.2|1,0><.1|2,0>}
    os << "{";
    for (MatingDistribution::Entries::const_iterator it=md.entries().begin(); it!=md.entries().end(); ++it)
        os << "<" << it->cumulativeWeight - (it==md.entries().begin()?0:(it-1)->cumulativeWeight)
           << "|" << it->indexPair.first << "," << it->indexPair.second << ">";
    os << "}";
    return os;
}


istream& operator>>(istream& is, MatingDistribution::Entry& entry)
{
    string buffer;
    getline(is, buffer, '>');
    if (!is) return is;

    char open, pipe, comma;

    istringstream iss(buffer);
    iss >> open >> entry.cumulativeWeight >> pipe >> entry.indexPair.first 
        >> comma >> entry.indexPair.second;

    if (!is || open!='<' || pipe!='|' || comma!=',')
        throw runtime_error("[operator>>(MatingDistribution::Entry)] Invalid input string.");

    return is;
}


istream& operator>>(istream& is, MatingDistribution& md)
{
    string buffer;
    getline(is, buffer,'}');
    if (!is) return is;

    char open;
    istringstream iss(buffer);
    iss >> open;
    if (open!='{')
        throw runtime_error("[operator>>(MatingDistribution)] Invalid input string.");

    // make a local copy of the entries, since cumulativeWeight isn't cumulative
    MatingDistribution::Entries entries;
    copy(istream_iterator<MatingDistribution::Entry>(iss), 
         istream_iterator<MatingDistribution::Entry>(),
         back_inserter(entries));     

    // clear any old data
    md = MatingDistribution();

    // md.push_back() does the cumulativeWeight calculation
    for (MatingDistribution::Entries::const_iterator it=entries.begin(); it!=entries.end(); ++it)
        md.push_back(it->cumulativeWeight, it->indexPair);

    return is;
}


//
// Population
//


namespace {


class RandomOrganismIndexGenerator
{
    public:

    RandomOrganismIndexGenerator(const Population& p,
                                 const DataVectorPtr& fitness_vector,
                                 const Random& random)
    :   population_size_(p.organisms().size()),
        fitness_cdf_max_(0),
        random_(random)
    {
        if (fitness_vector.get()) 
        {
            fitness_cdf_ = fitness_vector->cdf(); // memory allocation for cdf

            if (!fitness_cdf_.get() || fitness_cdf_->empty() || fitness_cdf_->size() != population_size_)
                throw runtime_error("[RandomOrganismIndexGenerator] This isn't happening.");

            fitness_cdf_max_ = fitness_cdf_->back();
        }
    }

    size_t operator()() const
    {
        if (!fitness_cdf_.get()) 
        {
            return random_.randint(0, population_size_-1); // uniform random index
        }
        else
        {
            // pick random index according to fitnesses
            double roll = random_.uniform(0, fitness_cdf_max_);
            DataVector::const_iterator it = lower_bound(fitness_cdf_->begin(), fitness_cdf_->end(), roll);
            return it - fitness_cdf_->begin();
        }
    }

    private:

    size_t population_size_;
    DataVectorPtr fitness_cdf_;
    double fitness_cdf_max_;
    const Random& random_;
};


typedef shared_ptr<RandomOrganismIndexGenerator> RandomOrganismIndexGeneratorPtr;


} // namespace


void Population::create_organisms(const Config& config,
                                  const Populations& populations,
                                  const DataVectorPtrs& fitnesses,
                                  const Random& random)
{
    if (config.size == 0)
        return;

    // sanity check: exactly one of {chromosomePairCount, matingDistribution} must be specified

    if (config.chromosomePairCount==0 && config.matingDistribution.empty())
        throw runtime_error("[Population::create_organisms()] Must specify nonzero chromosome pair count, or a mating distribution.");

    if (config.chromosomePairCount!=0 && !config.matingDistribution.empty())
        throw runtime_error("[Population::create_organisms()] Chromosome pair count and mating distribution both specified -- not sure what to do.");

    // create organisms from nothing

    if (config.chromosomePairCount != 0)
    {
        for (size_t i=0; i<config.size; ++i)
        {
            Chromosome::ID id(config.populationID, config.idOffset+i, 0, 0);
            organisms_.push_back(Organism(id, config.chromosomePairCount));
        }

        return;
    }

    // create organisms from previous generation

    if (populations.size() != fitnesses.size())
        throw runtime_error("[Population::create_organisms()] Fitness vector count != population count.");

    organisms_.reserve(config.size);

    // instantiate RandomOrganismIndexGenerators (one for each population)

    vector<RandomOrganismIndexGeneratorPtr> random_organism_index_generators;

    DataVectorPtrs::const_iterator fitness = fitnesses.begin();
    for (Populations::const_iterator population=populations.begin(); population!=populations.end(); ++population, ++fitness)
        random_organism_index_generators.push_back(RandomOrganismIndexGeneratorPtr(
            new RandomOrganismIndexGenerator(**population, *fitness, random)));

    // create Organisms for new population

    for (size_t i=0; i<config.size; ++i)
    {
        const MatingDistribution::IndexPair& parentIndices = config.matingDistribution.random_index_pair(random);

        if (max(parentIndices.first,parentIndices.second) >= populations.size())
            throw runtime_error("[Population::Population()] Indices out of bounds.");

        if (populations[parentIndices.first]->organisms().empty() ||
            populations[parentIndices.second]->organisms().empty())
            throw runtime_error("[Population::Population()] Empty population.");

        size_t index1 = (*random_organism_index_generators[parentIndices.first])();
        size_t index2 = 0;
        do { // avoid selfing
            index2 = (*random_organism_index_generators[parentIndices.second])();
        } while (parentIndices.first == parentIndices.second && index1 == index2);

        const Organism& mom = populations[parentIndices.first]->organisms()[index1];
        const Organism& dad = populations[parentIndices.second]->organisms()[index2];
        organisms_.push_back(Organism(mom, dad));
    }
}


shared_ptr<Population> Population::randomSubsample(size_t size, Random& random) const
{
    if (size > organisms_.size())
        throw runtime_error("[Population::randomSubsample] Sample size exceeds population size.");

    set<size_t> indices;
    while (indices.size() < size) // may take a long time if size is close to organisms_.size()
        indices.insert(random.randint(0,organisms_.size()-1));

    shared_ptr<Population> subsample(new Population());

    for (set<size_t>::const_iterator it=indices.begin(); it!=indices.end(); ++it)
        subsample->organisms_.push_back(organisms_[*it]);

    return subsample;
}


void Population::read(istream& is)
{
    size_t organism_count = 0;
    is.read((char*)&organism_count, sizeof(size_t));
    if (organism_count > 1e8) throw runtime_error("[Population::read()] Bad organism count.");
    organisms_.resize(organism_count);
    for (Organisms::iterator it=organisms_.begin(); it!=organisms_.end(); ++it)
        it->read(is);
}


void Population::write(ostream& os) const
{
    size_t organism_count = organisms_.size();
    os.write((const char*)&organism_count, sizeof(size_t));
    for (Organisms::const_iterator it=organisms_.begin(); it!=organisms_.end(); ++it)
        it->write(os);
}


PopulationsPtr Population::create_populations(const vector<Population::Config>& configs,
                                              const Populations& previous, 
                                              const DataVectorPtrs& fitnesses,
                                              const Random& random)
{
    PopulationsPtr result(new Populations);

    for (vector<Population::Config>::const_iterator it=configs.begin(); it!=configs.end(); ++it)
    {
        PopulationPtr p(new Population);
        p->create_organisms(*it, previous, fitnesses, random);
        result->push_back(p);
    }        

    return result;
}


bool operator==(const Population::Config& a, const Population::Config& b)
{
    return a.size == b.size &&
           a.populationID == b.populationID &&
           a.idOffset == b.idOffset &&
           a.chromosomePairCount == b.chromosomePairCount &&
           a.matingDistribution == b.matingDistribution;
}


bool operator!=(const Population::Config& a, const Population::Config& b)
{
    return !(a==b);
}


ostream& operator<<(ostream& os, const Population::Config& config)
{
    os << "population"
       << " size=" << config.size
       << " populationID=" << config.populationID;

    if (config.idOffset != 0)
        os << " idOffset=" << config.idOffset;

    if (config.chromosomePairCount != 0)
        os << " chromosomePairCount=" << config.chromosomePairCount;

    if (!config.matingDistribution.entries().empty())
        os << " matingDistribution=" << config.matingDistribution;

    return os;
}


istream& operator>>(istream& is, Population::Config& config)
{
    // population size=4 populationID=1 idOffset=1000 chromosomePairCount=3 matingDistribution={<0.42|0,0><0.66|1,0><0.23|2,1>}

    string buffer;
    getline(is, buffer);

    vector<string> tokens;
    istringstream iss(buffer);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

    if (tokens.empty() || tokens[0]!="population")
    {
        cerr << "--> " << buffer << endl;
        throw runtime_error("[operator<<(Population::Config)] Invalid config format");
    }

    for (vector<string>::const_iterator it=tokens.begin()+1; it!=tokens.end(); ++it)
    {
        size_t index_equal = it->find('=');
        if (index_equal == string::npos)
        {
            cerr << "Ignoring invalid token: " << *it << endl;
            continue;
        }

        string name = it->substr(0,index_equal);
        istringstream value(it->substr(index_equal+1));

        if (name.empty() || value.str().empty())
        {
            cerr << "Ignoring invalid token: " << *it << endl;
            continue;
        }

        if (name == "size")
            value >> config.size;
        else if (name == "populationID")
            value >> config.populationID;
        else if (name == "idOffset")
            value >> config.idOffset;
        else if (name == "chromosomePairCount")
            value >> config.chromosomePairCount;
        else if (name == "matingDistribution")
            value >> config.matingDistribution;
    }

    return is;
}


bool operator==(const Population& a, const Population& b)
{
    if (a.organisms().size() != b.organisms().size()) return false;

    for (vector<Organism>::const_iterator it=a.organisms().begin(), jt=b.organisms().begin(); 
         it!=a.organisms().end(); ++it, ++jt)
        if (*it != *jt) return false;

    return true;
}


bool operator!=(const Population& a, const Population& b)
{
    return !(a==b);
}


ostream& operator<<(ostream& os, const Population& p)
{
    for (vector<Organism>::const_iterator it=p.organisms().begin(); it!=p.organisms().end(); ++it)
        os << *it << endl;
    return os;
}


istream& operator>>(istream& is, Population& p)
{
    p.organisms_.clear();
    copy(istream_iterator<Organism>(is), istream_iterator<Organism>(), back_inserter(p.organisms_));
    return is;
}


