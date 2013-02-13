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


#ifndef _SIMULATIONCONTROLLER_NEUTRALADMIXTURE_HPP_
#define _SIMULATIONCONTROLLER_NEUTRALADMIXTURE_HPP_


#include "SimulationController.hpp"


class SimulationController_NeutralAdmixture : public SimulationController
{
    public:

    struct Config
    {
        std::string population_config_filename; // "popconfig"
        std::string output_directory;           // "outdir"

        Config(const Parameters& parameters = Parameters()); // allows auto conversion: Parameters->Config
    };

    SimulationController_NeutralAdmixture(const Config& config = Config());

    virtual void initialize();
    virtual void run() const;
    virtual void report() const;
    virtual void example(const std::string& output_directory) const;

    private:

    Config config_;
    Simulator::Config simulator_config_;
};


#endif //  _SIMULATIONCONTROLLER_NEUTRALADMIXTURE_HPP_

