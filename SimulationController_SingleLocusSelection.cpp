//
// SimulationController_SingleLocusSelection.cpp
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


#include "SimulationController_SingleLocusSelection.hpp"
#include "boost/filesystem.hpp"
#include <iostream>
#include <sstream>


using namespace std;
namespace bfs = boost::filesystem;


class Reporter_Temp : public Reporter // simple reporter for testing
{
    public:

    Reporter_Temp(const string& output_directory)
    :   outdir_(output_directory)
    {
        os_log_.open(outdir_ / "log.txt");
        if (!os_log_)
            throw runtime_error("[Reporter_Temp] Unable to open log.");
    }

    virtual void update(size_t generation_number,
                        const PopulationPtrs& populations,
                        const PopulationDatas& population_datas)
    {
        os_log_ << "generation: " << generation_number << endl;
        os_log_ << "populations: " << populations.size() << endl;
        assert(populations.size() == 1);
        os_log_ << "population 0 size: " << populations.front()->size() << endl;

        ostringstream filename;
        filename << "gen" << generation_number << ".txt"; 
        bfs::ofstream os_pop(outdir_ / filename.str());
        if (!os_pop)
            throw runtime_error(("[Reporter_Temp] Unable to open " + filename.str()).c_str());
        os_pop << *populations.front();
    }


    virtual void update_final(size_t generation_number,
                              const PopulationPtrs& populations,
                              const PopulationDatas& population_datas)
    {
    }

    private:

    bfs::path outdir_;
    bfs::ofstream os_log_;
};




SimulationController_SingleLocusSelection::Config::Config(const Parameters& parameters)
{
    seed = parameters.count("seed") ? atoi(parameters.at("seed").c_str()) : 0;
    output_directory = parameters.count("outdir") ? parameters.at("outdir") : "";

    population_size = parameters.count("popsize") ? atoi(parameters.at("popsize").c_str()) : 0;
    generation_count = parameters.count("gencount") ? atoi(parameters.at("gencount").c_str()) : 0;
    initial_allele_frequency = parameters.count("allelefreq") ? atoi(parameters.at("allelefreq").c_str()) : 0;

    w0 = parameters.count("w0") ? atof(parameters.at("w0").c_str()) : 1;
    w1 = parameters.count("w1") ? atof(parameters.at("w1").c_str()) : 1;
    w2 = parameters.count("w2") ? atof(parameters.at("w2").c_str()) : 1;
}


SimulationController_SingleLocusSelection::SimulationController_SingleLocusSelection(const Config& config)
:   config_(config)
{}


void SimulationController_SingleLocusSelection::SimulationController_SingleLocusSelection::usage() const
{   
    cout << "Usage:  simrecomb single_locus_selection [parameter_name=value] ...\n";
    cout << "        simrecomb sls [parameter_name=value] ...\n";
    cout << "        simrecomb sls config=config_filename ...\n";
    cout << endl;
    cout << "Required parameters:\n";
    cout << "  outdir=<output_directory>\n";
    cout << "  popsize=<population_size>\n";
    cout << "  gencount=<generation_count>\n";
    cout << endl;
    cout << "Optional parameters:\n";
    cout << "  config=<config_filename>\n";
    cout << "  seed=<value>\n";
    cout << "  allelefreq=<initial_allele_frequency> (default: allelefreq=0)\n";
    cout << "  w0=<relative_fitness_genotype_0> (default: w0=1)\n";
    cout << "  w1=<relative_fitness_genotype_1> (default: w1=1)\n";
    cout << "  w2=<relative_fitness_genotype_2> (default: w2=1)\n";
    cout << endl;
}


Random random_hack_; // TODO: remove


void SimulationController_SingleLocusSelection::initialize()
{
    // check parameters

    if (config_.output_directory.empty())
        throw runtime_error("[SimulationController_SingleLocusSelection] No output directory specified (outdir=value).");

    if (bfs::exists(config_.output_directory))
        throw runtime_error(("[SimulationController_SingleLocusSelection] Output directory exists: " + config_.output_directory).c_str());

    if (config_.population_size == 0)
        throw runtime_error("[SimulationController_SingleLocusSelection] Population size not specified (popsize=value).");

    if (config_.generation_count == 0)
        throw runtime_error("[SimulationController_SingleLocusSelection] Generation count not specified (gencount=value).");

    bfs::create_directories(config_.output_directory);

    // initialize simulator

    simulator_config_.seed = config_.seed;
    
    cout << "seed: " << config_.seed << endl; // TODO: remove

    simulator_config_.output_directory = config_.output_directory;    

    Population::Configs configs_gen_0(1);
    Population::Config& config_gen_0 = configs_gen_0.back();
    config_gen_0.size = config_.population_size;
    config_gen_0.populationID = 0; 
    config_gen_0.chromosomePairCount = 1; 
    simulator_config_.population_configs.push_back(configs_gen_0);

    Population::Configs configs_gen_next(1);
    Population::Config& config_gen_next = configs_gen_next.back();
    config_gen_next.size = config_.population_size;
    config_gen_next.populationID = 0; 
    config_gen_next.matingDistribution.push_back(1, make_pair(0,0));
    for (size_t i=0; i<config_.generation_count; ++i)
        simulator_config_.population_configs.push_back(configs_gen_next);

    cout << "ok\n";

    // SNPIndicator & popconfigs

    // QT_SingleLocusFitness
    // FF_Identity

    // Reporters:  
    //   full population for debugging 
    //   allele freq / homozygosity
    //   block lengths?
    //   mean fitness

    simulator_config_.reporters.push_back(ReporterPtr(new Reporter_Temp(config_.output_directory)));

    simulator_ = SimulatorPtr(new Simulator(simulator_config_));
}


void SimulationController_SingleLocusSelection::run() const
{
    // hack: TODO remove -- fix Simulator logic for this

    Organism::recombinationPositionGenerator_ =
        shared_ptr<RecombinationPositionGenerator>(new RecombinationPositionGenerator_Trivial(random_hack_));

    simulator_->simulate_all();
}


void SimulationController_SingleLocusSelection::report() const
{
    simulator_->update_final();
}


void SimulationController_SingleLocusSelection::example(const std::string& output_directory) const
{
    cout << "SimulationController_SingleLocusSelection::example() not implemented.\n"; // TODO
}



