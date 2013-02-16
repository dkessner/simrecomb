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
#include <sstream>
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"


using namespace std;
namespace bfs = boost::filesystem;


//
// Simulator
//


Simulator::Simulator(const Config& config)
:   config_(config), random_(config.seed),
    current_generation_(0), current_populations_(new Populations)
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


void Simulator::simulate_single_generation(ostream* os_log)
{
    if (!current_populations_.get())
        throw runtime_error("[Simulator::simulate_single_generation()] Null pointer.");

    if (current_generation_ >= config_.population_configs.size())
        throw runtime_error("[Simulator::simulate_single_generation()] Population config not specified.");

    if (config_.os_progress) *config_.os_progress << "[Simulator] Generation " << current_generation_ << endl;

    DataVectorPtrs dummy_fitnesses(current_populations_->size()); 

    PopulationsPtr next = Population::create_populations(
        config_.population_configs[current_generation_], *current_populations_, dummy_fitnesses, random_);

    if (os_log) *os_log << current_generation_ << endl;

    current_populations_ = next;
    ++current_generation_;
}


void Simulator::simulate_all()
{
    bfs::path output_directory(config_.output_directory);

    bfs::ofstream os_log(output_directory / "log.txt");

    const size_t generation_count = config_.population_configs.size();
    for (size_t generation=0; generation<generation_count; generation++)
    {
        simulate_single_generation(&os_log);
    }

    os_log.close();

    // output the last generation
    
    for (size_t i=0; i<current_populations_->size(); ++i)
    {
        ostringstream filename;
        filename << "pop" << i << ".txt"; 
        bfs::ofstream os_pop(output_directory / filename.str());
        os_pop << *(*current_populations_)[i];
    }
}


