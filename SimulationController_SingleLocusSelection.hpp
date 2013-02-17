//
// SimulationController_SingleLocusSelection.hpp
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


#ifndef _SIMULATIONCONTROLLER_SINGLELOCUSSELECTION_HPP_
#define _SIMULATIONCONTROLLER_SINGLELOCUSSELECTION_HPP_


#include "SimulationController.hpp"


class SimulationController_SingleLocusSelection : public SimulationController
{
    public:

    struct Config
    {
        int seed;
        std::string output_directory;           // "outdir"

        size_t population_size;                 // "popsize"
        size_t generation_count;                // "gencount"
        size_t initial_allele_frequency;        // "allelefreq"
        
        // relative fitnesses for genotype g in {0,1,2}, where g is the # of selected alleles
        std::vector<double> w;                       // "w0", "w1", "w2"

        Config(const Parameters& parameters = Parameters()); // allows auto conversion: Parameters->Config
    };

    SimulationController_SingleLocusSelection(const Config& config = Config());

    virtual void usage() const;
    virtual void initialize();
    virtual void run() const;
    virtual void report() const;
    virtual void example(const std::string& output_directory) const;

    private:

    Config config_;
    Simulator::Config simulator_config_;
    SimulatorPtr simulator_;
};


#endif //  _SIMULATIONCONTROLLER_SINGLELOCUSSELECTION_HPP_

