//
// PopulationAnalyzer.hpp
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


#ifndef _POPULATIONANALYZER_HPP_
#define _POPULATIONANALYZER_HPP_


#include <vector>
#include <map>
#include "boost/shared_ptr.hpp"
#include "Population.hpp"


using boost::shared_ptr;


class DataVector : public std::vector<double>
{
    // add convenience functions:  sum, sum_squares, mean, variance, CDF
};


typedef shared_ptr<DataVector> DataVectorPtr;


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


struct PopulationData
{
    std::map<Locus, DataVectorPtr> genotypes;
    std::vector<DataVectorPtr> trait_values;
    DataVectorPtr fitness;
};


class PopulationAnalyzer
{
    public:

    shared_ptr<PopulationData> analyze(const Population& population) const;
    
    private:

    // list of QuantitativeTraits
    // ref to Genotyper
    // ref to FitnessFunction
    // configuration flags for output, and output functions (or delegate)
};


#endif //  _POPULATIONANALYZER_HPP_

