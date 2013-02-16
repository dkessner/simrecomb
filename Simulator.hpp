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


#include "QuantitativeTrait.hpp"
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
        unsigned int seed;                                      // for Random
        std::string output_directory;                           // all output files placed here
        std::ostream* os_progress;                              // progress update stream (default: stdout)

        std::vector<std::string> genetic_map_filenames;         // one filename for each chromosome pair
        std::vector<Population::Configs> population_configs;    // Population::Configs for each generation

        SNPIndicatorPtr snp_indicator;
        QuantitativeTraitPtrs quantitative_traits;
        FitnessFunctionPtr fitness_function;
        ReporterPtrs reporters;

        Config() : seed(0), os_progress(&std::cout) {}
    };

    Simulator(const Config& config);  

    void simulate_single_generation(std::ostream* os_log = 0);
    void simulate_all();

    private:

    Config config_;
    Random random_;
    Genotyper genotyper_;

    size_t current_generation_;
    PopulationsPtr current_populations_;
    PopulationDatasPtr current_population_datas_;
};


typedef shared_ptr<Simulator> SimulatorPtr;


#endif //  _SIMULATOR_HPP_

