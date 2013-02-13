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
    population_config_filename = parameters.at("popconfig");
    output_directory = parameters.at("outdir");
}


SimulationController_NeutralAdmixture::SimulationController_NeutralAdmixture(const Config& config)
:   config_(config)
{
    cout << "we're here\n";
}


void SimulationController_NeutralAdmixture::example() const
{
}


void SimulationController_NeutralAdmixture::initialize()
{
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




