//
// SimulationController_NeutralAdmixture.hpp
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


#include "SimulationController_NeutralAdmixture.hpp"
#include "boost/filesystem.hpp"
#include <iostream>


using namespace std;
namespace bfs = boost::filesystem;


SimulationController_NeutralAdmixture::Config::Config(const Parameters& parameters)
{
    population_config_filename = parameters.count("popconfig") ? parameters.at("popconfig") : "";
    output_directory = parameters.count("outdir") ? parameters.at("outdir") : "";
}


SimulationController_NeutralAdmixture::SimulationController_NeutralAdmixture(const Config& config)
:   config_(config)
{}


void SimulationController_NeutralAdmixture::example(const string& output_directory) const
{   
    Simulator::Config config;

    const size_t chromosomePairCount_ = 3;
    const unsigned int populationSize_ = 10000;
    const double admixtureProportion_ = .8; // fraction of genes from 1st population

    for (size_t i=0; i<chromosomePairCount_; i++) 
        config.geneticMapFilenames.push_back("genetic_map_chr21_b36.txt"); // hack

    // generation 0 (ancestral populations)

    config.populationConfigs.push_back(vector<Population::Config>(3));

    Population::Config* config_pop = &config.populationConfigs[0][0];
    config_pop->size = 0;
    config_pop->populationID = 0;

    config_pop = &config.populationConfigs[0][1];
    config_pop->size = populationSize_;
    config_pop->populationID = 1;
    config_pop->chromosomePairCount = chromosomePairCount_;

    config_pop = &config.populationConfigs[0][2];
    config_pop->size = populationSize_;
    config_pop->populationID = 2;
    config_pop->chromosomePairCount = chromosomePairCount_;

    // generation 1 (initial admixture)

    config.populationConfigs.push_back(vector<Population::Config>(1));

    config_pop = &config.populationConfigs[1][0];
    config_pop->size = populationSize_;
    double p = admixtureProportion_;
    config_pop->matingDistribution.push_back(p*p, make_pair(1,1));
    config_pop->matingDistribution.push_back(2*p*(1-p), make_pair(1,2));
    config_pop->matingDistribution.push_back((1-p)*(1-p), make_pair(2,2));

    // subsequent generations - just recombination

    for (size_t generation=2; generation<8; generation++)
    {
        config.populationConfigs.push_back(vector<Population::Config>(1));
        config_pop = &config.populationConfigs[generation][0];
        config_pop->size = populationSize_;
        config_pop->matingDistribution.push_back(1, make_pair(0,0));
    }

    cout << config;
}


void SimulationController_NeutralAdmixture::initialize()
{
    if (config_.output_directory.empty())
        throw runtime_error("[SimulationController_NeutralAdmixture] No output directory specified (outdir=value).");

    if (config_.population_config_filename.empty())
        throw runtime_error("[SimulationController_NeutralAdmixture] No population config file specified (popconfig=value).");

    if (!bfs::exists(config_.population_config_filename))
        throw runtime_error(("[SimulationController_NeutralAdmixture] Population config file not found: " + config_.population_config_filename).c_str());

    if (bfs::exists(config_.output_directory))
        throw runtime_error(("[SimulationController_NeutralAdmixture] Output directory exists: " + config_.output_directory).c_str());

    cout << "[SimulationController_NeutralAdmixture] Reading population configuration file " << config_.population_config_filename << endl;
    bfs::ifstream is(config_.population_config_filename);
    is >> simulator_config_;
    is.close();

    bfs::create_directories(config_.output_directory);
}


void SimulationController_NeutralAdmixture::run() const
{
    Simulator simulator(simulator_config_, config_.output_directory);
    simulator.simulate_all();
}


void SimulationController_NeutralAdmixture::report() const
{
}


