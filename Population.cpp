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


Population::Population(const Config& config, const Random& random)
:   random_(random)
{
    for (size_t i=0; i<config.size; ++i)
    {
        Chromosome::ID id(config.populationID, config.idOffset+i, 0, 0);
        organisms_.push_back(Organism(id, config.chromosomePairCount));
    }
}


Population::Population(const Config& config,
                       const Populations& populations,
                       const Random& random)
:   random_(random)
{
    for (size_t i=0; i<config.size; ++i)
    {
        const MatingDistribution::IndexPair& parentIndices = config.matingDistribution.random_index_pair(random_);

        if (max(parentIndices.first,parentIndices.second) >= populations.size())
            throw runtime_error("[Population::Population()] Indices out of bounds.");

        if (populations[parentIndices.first]->organisms().empty() ||
            populations[parentIndices.second]->organisms().empty())
            throw runtime_error("[Population::Population()] Empty population.");

        size_t index1 = random_.randint(0, populations[parentIndices.first]->organisms().size()-1);
        size_t index2 = 0;
        do { // avoid selfing
            index2 = random_.randint(0, populations[parentIndices.second]->organisms().size()-1);
        } while (parentIndices.first == parentIndices.second && index1 == index2);

        const Organism& mom = populations[parentIndices.first]->organisms()[index1];
        const Organism& dad = populations[parentIndices.second]->organisms()[index2];

        Organism::Gamete egg = mom.create_gamete();
        Organism::Gamete sperm = dad.create_gamete();
        organisms_.push_back(Organism(egg, sperm));
    }
}


Population::Population(const std::string& filename, const Random& random)
:   random_(random)
{
    ifstream is(filename.c_str());
    if (!is)
        throw runtime_error(("[Population] Unable to open file " + filename).c_str());
        
    is >> *this;

    if (organisms_.empty())
        cerr << "[Population] Warning: no population data read from file " << filename << endl;
}


shared_ptr<Population> Population::randomSubsample(size_t size) const
{
    if (size > organisms_.size())
        throw runtime_error("[Population::randomSubsample] Sample size exceeds population size.");

    set<size_t> indices;
    while (indices.size() < size) // may take a long time if size is close to organisms_.size()
        indices.insert(random_.randint(0,organisms_.size()-1));

    shared_ptr<Population> subsample(new Population(Config(), random_));

    for (set<size_t>::const_iterator it=indices.begin(); it!=indices.end(); ++it)
        subsample->organisms_.push_back(organisms_[*it]);

    return subsample;
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

    if (config.chromosomePairCount != 1)
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


