//
// Simulator.hpp
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


#ifndef _SIMULATOR_HPP_
#define _SIMULATOR_HPP_


#include "Population.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include <vector>
#include <string>
#include <iostream>


class Simulator
{
    public:

    struct Config
    {
        unsigned int seed;
        std::vector<std::string> geneticMapFilenames;     // one filename for each chromosome pair
        std::vector<Population::Configs> populationConfigs; // Population::Configs for each generation

        Config() : seed(0) {}
    };

    Simulator(const Config& config, 
              const std::string& output_directory,      // all output placed here
              std::ostream* os_progress = &std::cout);  // progress update stream (default: stdout)

    void simulate_single_generation(std::ostream* os_log = 0);
    void simulate_all();

    private:

    Config config_;
    Random random_;

    boost::filesystem::path output_directory_;
    std::ostream* os_progress_;

    size_t current_generation_;
    PopulationsPtr current_populations_;
};


#endif //  _SIMULATOR_HPP_

