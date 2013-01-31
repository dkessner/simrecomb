//
// Genotyper.hpp
//
// Copyright 2013 Darren Kessner
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


#ifndef _GENOTYPER_HPP_
#define _GENOTYPER_HPP_


#include "Organism.hpp"


struct Locus
{
    size_t chromosome_pair_index;
    unsigned int position;

    Locus(size_t _chromosome_pair_index = 0, unsigned int _position = 0) 
    :   chromosome_pair_index(_chromosome_pair_index), position(_position) 
    {}
};


inline bool operator<(const Locus& a, const Locus& b)
{
    return (a.chromosome_pair_index < b.chromosome_pair_index) ||
           (a.chromosome_pair_index == b.chromosome_pair_index && a.position < b.position);
}


/*
struct PopulationData
{
    std::map<Locus, DataVectorPtr> genotypes;
    std::vector<DataVectorPtr> trait_values;
    DataVectorPtr fitness;
};
*/


//////////////



class SNPIndicator
{
    public:

    // return value in {0, 1}
    virtual unsigned int operator()(unsigned int position, unsigned int chromosome_id) const = 0;
    virtual ~SNPIndicator() {}
};


class Genotyper
{
    public:

    // returns genotype value in {0, 1, 2}
    unsigned int genotype(const Organism& organism, size_t chromosome_pair_index,
                          unsigned int position, const SNPIndicator& indicator) const;

    // returns (population x position) genotype matrix, filled in via iteration of genotype(organism, position)
    // multi_array genotype(const Population& population, vector<unsigned int> position, 
    //                      const SNPIndicator& indicator) const;
};


#endif //  _GENOTYPER_HPP_

