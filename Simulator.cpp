//
// Simulator.cpp
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


#include "Simulator.hpp"
#include <iostream>
#include <iterator>
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"


using namespace std;
namespace bfs = boost::filesystem;


//
// Simulator
//


Simulator::Simulator(const Config& config)
:   config_(config), 
    random_(config.seed),
    current_generation_(0), 
    current_populations_(new PopulationPtrs),
    current_population_datas_(new PopulationDatas)
{
    cout << "[Simulator] Initializing.\n";

    // use trivial implementations by default

    if (!config_.snp_indicator.get())
        config_.snp_indicator = SNPIndicatorPtr(new SNPIndicator_Trivial);

    if (!config_.fitness_function.get())
        config_.fitness_function = FitnessFunctionPtr(new FitnessFunction_Trivial);

    // initialize recombination maps

    cout << "[Simulator] Initializing recombination maps.\n";
    Organism::recombinationPositionGenerator_ =
        shared_ptr<RecombinationPositionGenerator>(
            new RecombinationPositionGenerator_RecombinationMap(config.genetic_map_filenames, random_));
}


void Simulator::simulate_single_generation()
{
    // sanity checks

    if (!current_populations_.get())
        throw runtime_error("[Simulator::simulate_single_generation()] Null pointer.");

    if (current_generation_ >= config_.population_configs.size())
        throw runtime_error("[Simulator::simulate_single_generation()] Population config not specified.");

    if (config_.os_progress) *config_.os_progress << "[Simulator] Generation " << current_generation_ << endl;

    // create next generation

    DataVectorPtrs fitnesses;
    for (PopulationDatas::const_iterator popdata=current_population_datas_->begin();
         popdata!=current_population_datas_->end(); ++popdata)
        fitnesses.push_back(popdata->fitnesses);

    PopulationPtrsPtr next_populations = Population::create_populations(
        config_.population_configs[current_generation_], *current_populations_, fitnesses, random_);

    // collect data on the populations

    PopulationDatasPtr next_population_datas(new PopulationDatas(next_populations->size()));

    // calculate genotypes

    Loci loci_all;

    for (QuantitativeTraitPtrs::const_iterator qt=config_.quantitative_traits.begin();
         qt!=config_.quantitative_traits.end(); ++qt)
    {
        const Loci& loci = (*qt)->loci();
        for (Loci::const_iterator locus=loci.begin(); locus!=loci.end(); ++locus)
            loci_all.insert(*locus);
    }

    PopulationPtrs::const_iterator population = next_populations->begin();
    for (PopulationDatas::iterator popdata=next_population_datas->begin();
         popdata!=next_population_datas->end(); ++popdata, ++population)
    {
         popdata->genotypes = genotyper_.genotype(loci_all, **population, *config_.snp_indicator);
    }

    // calculate quantitative trait values

    for (PopulationDatas::iterator popdata=next_population_datas->begin();
         popdata!=next_population_datas->end(); ++popdata)
    {
        for (QuantitativeTraitPtrs::const_iterator qt=config_.quantitative_traits.begin();
             qt!=config_.quantitative_traits.end(); ++qt)
        {
            (*popdata->trait_values)[(*qt)->id()] = (*qt)->calculate_trait_values(popdata->genotypes);            
        }
    }

    // calculate fitnesses

    for (PopulationDatas::iterator popdata=next_population_datas->begin();
         popdata!=next_population_datas->end(); ++popdata)
    {
         if (popdata->trait_values.get())
             popdata->fitnesses = config_.fitness_function->calculate_fitness(*popdata->trait_values);
    }

    // update

    current_populations_ = next_populations;
    current_population_datas_ = next_population_datas;

    for (ReporterPtrs::iterator reporter=config_.reporters.begin(); reporter!=config_.reporters.end(); ++reporter)
    {
        (*reporter)->update(current_generation_, *current_populations_, *current_population_datas_);
    }

    ++current_generation_;
}


void Simulator::simulate_all()
{
    const size_t generation_count = config_.population_configs.size();
    for (size_t generation=0; generation<generation_count; generation++)
    {
        simulate_single_generation();
    }
}


void Simulator::update_final()
{
    for (ReporterPtrs::iterator reporter=config_.reporters.begin(); reporter!=config_.reporters.end(); ++reporter)
    {
        (*reporter)->update_final(current_generation_, *current_populations_, *current_population_datas_);
    }
}


