//
// QuantitativeTrait.hpp
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


#ifndef _QUANTITATIVETRAIT_HPP_
#define _QUANTITATIVETRAIT_HPP_


#include "DataVector.hpp"


class QuantitativeTrait
{
    public:

    // return list of loci -- used by Genotyper
    virtual const Loci& loci() const = 0;

    // calculate trait values for a single population using genotype map
    virtual DataVectorPtr calculate_trait_values(GenotypeMapPtr genotype_map) const = 0;

    virtual ~QuantitativeTrait() {}
};


class QuantitativeTrait_SingleLocusFitness : public QuantitativeTrait
{
    public:

    virtual const Loci& loci() const {return loci_;}

    // transform {0, 1, 2} -> {1, 1+hs, 1+2s}
    virtual DataVectorPtr calculate_trait_values(GenotypeMapPtr genotype_map) const;

    private:

    Loci loci_; // vector with a single locus
    double s_; // selection coefficient
    double h_; // dominance
};


class FitnessFunction
{
    public:

    DataVectorPtr calculate_fitness(const DataVectorPtrs& trait_values) const = 0;

    virtual ~FitnessFunction() {}
};


class FitnessFunction_Identity : public FitnessFunction
{
    public:

    DataVectorPtr calculate_fitness(const DataVectorPtrs& trait_values) const
    {
        assert(trait_values.size() == 1);  // -> runtime_error
        return trait_values[0];
    }
};


#endif //  _QUANTITATIVETRAIT_HPP_

