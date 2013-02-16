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


#include "DataVector.hpp"
#include "Organism.hpp"
#include "shared_ptr.hpp"
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
    bool empty() const {return entries_.empty();} 

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
typedef shared_ptr<Population> PopulationPtr;
typedef std::vector<PopulationPtr> PopulationPtrs;
typedef shared_ptr<PopulationPtrs> PopulationPtrsPtr;


class Population
{
    public:

    Population() {}

    Population(const Organisms& organisms)
    :   organisms_(organisms)
    {}

    struct Config
    {
        size_t size;

        // for generating population from nothing
        size_t chromosomePairCount; 
        unsigned int populationID;
        unsigned int idOffset;

        // for generating population from a previous generation
        MatingDistribution matingDistribution;

        Config() 
        :   size(0), chromosomePairCount(0), populationID(0), idOffset(0)
        {}
    };

    typedef std::vector<Config> Configs;

    void create_organisms(const Config& config,
                          const PopulationPtrs& populations = PopulationPtrs(),
                          const DataVectorPtrs& fitnesses = DataVectorPtrs(), // null ok, but size must match populations
                          const Random& random = Random());

    const std::vector<Organism>& organisms() const {return organisms_;}
    size_t size() const {return organisms_.size();}

    shared_ptr<Population> randomSubsample(size_t size, Random& random) const;

    // binary read/write

    void read(std::istream& is);
    void write(std::ostream& os) const;

    // convenience function: creates new generation from previous

    static PopulationPtrsPtr create_populations(const std::vector<Population::Config>& configs,
                                                const PopulationPtrs& previous, 
                                                const DataVectorPtrs& fitnesses,
                                                const Random& random);

    private:

    Organisms organisms_;

    friend std::istream& operator>>(std::istream& is, Population& p);

    // disallow copying
    Population(Population&);
    Population& operator=(Population&);
};


bool operator==(const Population::Config& a, const Population::Config& b);
bool operator!=(const Population::Config& a, const Population::Config& b);
std::ostream& operator<<(std::ostream& os, const Population::Config& config);
std::istream& operator>>(std::istream& is, Population::Config& config);


bool operator==(const Population& a, const Population& b);
bool operator!=(const Population& a, const Population& b);
std::ostream& operator<<(std::ostream& os, const Population& p);
std::istream& operator>>(std::istream& is, Population& p);


// multiple generations: each generation has Population::Configs, one Config per population
std::ostream& operator<<(std::ostream& os, const std::vector<Population::Configs>& generation_configs);
std::istream& operator>>(std::istream& is, std::vector<Population::Configs>& generation_configs);


#endif // _POPULATION_HPP_

