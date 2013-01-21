//
// Organism.cpp
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

#include "Organism.hpp"
#include "Random.hpp"
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <iterator>


using namespace std;
using boost::shared_ptr;


//
// RecombinationPositionGenerator
//

vector<unsigned int> RecombinationPositionGenerator_Trivial::get_positions(size_t index) const
{
    vector<unsigned int> result;
    if (random_.random()>=.5) result.push_back(0); // start with 2nd chromosome
    return result;
}


vector<unsigned int> RecombinationPositionGenerator_RecombinationMap::get_positions(size_t index) const
{
    if (index >= recombinationMaps_.size())
        throw runtime_error("[RecombinationPositionGenerator_RecombinationMap::get_positions()] Index out of bounds.");

    vector<unsigned int> positions = recombinationMaps_[index]->random_positions(); // may not be sorted 
    if (random_.random()>=.5) positions.push_back(0); // start with 2nd chromosome
    sort(positions.begin(), positions.end());

    return positions; 
}


RecombinationPositionGenerator_RecombinationMap::RecombinationPositionGenerator_RecombinationMap(
        const vector<string>& filenames,
        const Random& random)
:   random_(random)
{
    for (vector<string>::const_iterator it=filenames.begin(); it!=filenames.end(); ++it)
        recombinationMaps_.push_back(shared_ptr<RecombinationMap>(
            new RecombinationMap(*it, random)));
}


//
// Organism
//


shared_ptr<RecombinationPositionGenerator> Organism::recombinationPositionGenerator_; // static storage


Organism::Organism(unsigned int id, size_t chromosomeCount)
{
    Chromosome::ID id0(id);
    Chromosome::ID id1(id);
    id0.which = 0;
    id1.which = 1;

    for (size_t i=0; i<chromosomeCount; ++i)
    {
        id0.pair = id1.pair = i;
        chromosomePairs_.push_back(make_pair(Chromosome(id0), Chromosome(id1)));
    }
}


Organism::Organism(const Gamete& g1, const Gamete& g2)
{
    if (g1.size() != g2.size())
        throw runtime_error("[Organism::Organism()] Different gamete sizes.");

    for (size_t i=0; i<g1.size(); ++i)
        chromosomePairs_.push_back(make_pair(g1[i],g2[i]));
}


Organism::Organism(const Organism& mom, const Organism& dad)
{
    if (!recombinationPositionGenerator_.get())
        throw runtime_error("[Organism::Organism(mom, dad)] No RecombinationPositionGenerator.");

    if (mom.chromosomePairs_.size() != dad.chromosomePairs_.size())
        throw runtime_error("[Organism::Organism(mom, dad)] Parents chromosome counts differ.");

    this->chromosomePairs_.reserve(mom.chromosomePairs_.size());

    size_t chromosome_index = 0;
    for (ChromosomePairs::const_iterator it=mom.chromosomePairs_.begin(), jt=dad.chromosomePairs_.begin();
         it!=mom.chromosomePairs_.end(); ++it, ++jt, ++chromosome_index)
    {
        vector<unsigned int> positions_mom = recombinationPositionGenerator_->get_positions(chromosome_index);
        vector<unsigned int> positions_dad = recombinationPositionGenerator_->get_positions(chromosome_index);

        this->chromosomePairs_.push_back(make_pair(
            Chromosome(it->first, it->second, positions_mom),
            Chromosome(jt->first, jt->second, positions_dad)));
    }
}


Organism::Gamete Organism::create_gamete() const
{
    if (!recombinationPositionGenerator_.get())
        throw runtime_error("[Organism::create_gamete()] No RecombinationPositionGenerator");

    Gamete result;

    for (ChromosomePairs::const_iterator it=chromosomePairs_.begin(); it!=chromosomePairs_.end(); ++it)
    {
        vector<unsigned int> positions = 
            recombinationPositionGenerator_->get_positions(it-chromosomePairs_.begin());

        result.push_back(Chromosome(it->first, it->second, positions));
    }

    return result;
}


bool operator==(const Organism& a, const Organism& b)
{
    if (a.chromosomePairs().size() != b.chromosomePairs().size())
        return false;

    for (ChromosomePairs::const_iterator it=a.chromosomePairs().begin(), jt=b.chromosomePairs().begin();
         it!=a.chromosomePairs().end(); ++it, ++jt)
        if (*it!=*jt) return false;

    return true;
}


bool operator!=(const Organism& a, const Organism& b)
{
    return !(a==b);
}


ostream& operator<<(ostream& os, const Organism& o)
{
    for (ChromosomePairs::const_iterator it=o.chromosomePairs().begin(); it!=o.chromosomePairs().end(); ++it)
        os << "+ " << it->first << "\n- " << it->second << endl;
    return os;
}


istream& operator>>(istream& is, Organism& o)
{
    Organism::Gamete g_plus, g_minus;

    while (is)
    {
        string buffer;
        getline(is, buffer);
        if (!is || buffer.empty()) break;

        char plus_minus;
        Chromosome chromosome(0);
        istringstream iss(buffer);
        iss >> plus_minus >> chromosome;

        if (plus_minus == '+')
            g_plus.push_back(chromosome); 
        else if (plus_minus == '-')
            g_minus.push_back(chromosome); 
        else
            throw runtime_error("[operator>>(Organism)] Invalid format.");
    }

    if (!g_plus.empty())
        o = Organism(g_plus, g_minus);

    return is;
}


