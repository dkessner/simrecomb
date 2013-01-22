//
// Population.hpp
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

#ifndef _POPULATION_HPP_
#define _POPULATION_HPP_


#include "Organism.hpp"
#include "boost/shared_ptr.hpp"
#include <vector>


class MatingDistribution
{
    public:

    typedef std::pair<size_t,size_t> IndexPair;

    MatingDistribution() : totalWeight_(0) {}
    void push_back(double weight, const IndexPair& indexPair);
    const IndexPair& random_index_pair(const Random& random) const;

    struct Entry
    {
        double cumulativeWeight;
        IndexPair indexPair;
        Entry(double c=0, const IndexPair& i = std::make_pair(0,0)) : cumulativeWeight(c), indexPair(i) {}
    };

    typedef std::vector<Entry> Entries;
    const Entries& entries() const {return entries_;} 

    private:
    Entries entries_;
    double totalWeight_;
};


bool operator==(const MatingDistribution::Entry& a, const MatingDistribution::Entry& b);
bool operator!=(const MatingDistribution::Entry& a, const MatingDistribution::Entry& b);
bool operator==(const MatingDistribution& a, const MatingDistribution& b);
bool operator!=(const MatingDistribution& a, const MatingDistribution& b);
std::ostream& operator<<(std::ostream& os, const MatingDistribution& md);
std::istream& operator>>(std::istream& is, MatingDistribution::Entry& entry);
std::istream& operator>>(std::istream& is, MatingDistribution& md);


class Population;
typedef std::vector< boost::shared_ptr<Population> > Populations;


class Population
{
    public:

    struct Config
    {
        size_t size;

        // for generating population from nothing
        unsigned int populationID;
        unsigned int idOffset;
        size_t chromosomePairCount; 

        // for generating population from a previous generation
        MatingDistribution matingDistribution;

        Config() 
        :   size(0), populationID(0), idOffset(0), chromosomePairCount(0)
        {}
    };

    // construct an initial Population
    Population(const Config& config = Config());

    // construct a Population from Populations (e.g. from a previous generation)
    Population(const Config& config,
               const Populations& populations,
               const Random& random);

    // construct from file
    //Population(const std::string& filename); // TODO: make binary or remove?

    const std::vector<Organism>& organisms() const {return organisms_;}

    boost::shared_ptr<Population> randomSubsample(size_t size, Random& random) const;

    // binary read/write
    void read(std::istream& is);
    void write(std::ostream& os) const;

    private:

    std::vector<Organism> organisms_;

    friend std::istream& operator>>(std::istream& is, Population& p);
};


bool operator==(const Population::Config& a, const Population::Config& b);
bool operator!=(const Population::Config& a, const Population::Config& b);
std::ostream& operator<<(std::ostream& os, const Population::Config& config);
std::istream& operator>>(std::istream& is, Population::Config& config);


bool operator==(const Population& a, const Population& b);
bool operator!=(const Population& a, const Population& b);
std::ostream& operator<<(std::ostream& os, const Population& p);
std::istream& operator>>(std::istream& is, Population& p);


#endif // _POPULATION_HPP_

