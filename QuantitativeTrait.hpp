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


#include "Genotyper.hpp"
#include "DataVector.hpp"
#include "shared_ptr.hpp"
#include <stdexcept>


class QuantitativeTrait
{
    public:

    QuantitativeTrait(int id)
    :   id_(id)
    {}

    // identifier for the trait -- should be set by whoever instantiates the QuantitativeTrait
    int id() const {return id_;}

    // return set of loci (the QTLs contributing to the trait)
    const Loci& loci() const {return loci_;}

    // calculate trait values for a single population using genotypes
    virtual DataVectorPtr calculate_trait_values(GenotypeMapPtr genotypes) const = 0;

    virtual ~QuantitativeTrait() {}

    protected:

    int id_;
    Loci loci_;
};


typedef shared_ptr<QuantitativeTrait> QuantitativeTraitPtr;
typedef std::vector<QuantitativeTraitPtr> QuantitativeTraitPtrs;


class QuantitativeTrait_SingleLocusFitness : public QuantitativeTrait
{
    public:

    QuantitativeTrait_SingleLocusFitness(int id)
    :   QuantitativeTrait(id)
    {}

    // transform {0, 1, 2} -> {1, 1+hs, 1+2s}
    virtual DataVectorPtr calculate_trait_values(GenotypeMapPtr genotype_map) const;

    private:

    double s_; // selection coefficient
    double h_; // dominance
};


typedef std::map<int, DataVectorPtr> TraitValueMap;  // map QT id -> trait_values
typedef shared_ptr<TraitValueMap> TraitValueMapPtr;


class FitnessFunction
{
    public:

    virtual DataVectorPtr calculate_fitness(const TraitValueMap& trait_values) const = 0;
    virtual ~FitnessFunction() {}
};


typedef shared_ptr<FitnessFunction> FitnessFunctionPtr;


class FitnessFunction_Trivial : public FitnessFunction
{
    public:

    virtual DataVectorPtr calculate_fitness(const TraitValueMap& trait_values) const
    {
        return DataVectorPtr();
    }
};


class FitnessFunction_Identity : public FitnessFunction
{
    public:

    FitnessFunction_Identity(int qtid)
    :   qtid_(qtid)
    {}

    DataVectorPtr calculate_fitness(const TraitValueMap& trait_values) const
    {
        if (!trait_values.count(qtid_))
            throw std::runtime_error("[FitnessFunction_Identity] Quantitative trait id not found.");
        return trait_values.at(qtid_);
    }

    private:

    int qtid_;
};


struct PopulationData
{
    GenotypeMapPtr genotypes;
    TraitValueMapPtr trait_values;
    DataVectorPtr fitnesses;
};


typedef std::vector<PopulationData> PopulationDatas;
typedef shared_ptr<PopulationDatas> PopulationDatasPtr;


class Reporter
{
    public:

    virtual void update(size_t generation_number,
                        const Populations& populations,
                        const PopulationDatas& population_datas) = 0;

    virtual ~Reporter(){}
};


typedef shared_ptr<Reporter> ReporterPtr;
typedef std::vector<ReporterPtr> ReporterPtrs;


#endif //  _QUANTITATIVETRAIT_HPP_

