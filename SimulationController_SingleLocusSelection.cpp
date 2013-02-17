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
#include "boost/lambda/lambda.hpp"
#include <iostream>
#include <sstream>


using namespace std;
namespace bfs = boost::filesystem;
using namespace boost::lambda;


//////////////////////////////////////////////////////////////////////////////////////////////////
//
// TODO:  move this stuff

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
        filename << "pop0_" << generation_number << ".txt"; 
        bfs::ofstream os_pop(outdir_ / filename.str());
        if (!os_pop)
            throw runtime_error(("[Reporter_Temp] Unable to open " + filename.str()).c_str());
        os_pop << *populations.front();
    }

    virtual void update_final(size_t generation_number,
                              const PopulationPtrs& populations,
                              const PopulationDatas& population_datas)
    {}

    private:

    bfs::path outdir_;
    bfs::ofstream os_log_;
};


class Reporter_Genotypes : public Reporter
{
    public:

    Reporter_Genotypes(const string& output_directory, Locus locus)
    :   outdir_(output_directory), locus_(locus)
    {}

    virtual void update(size_t generation_number,
                        const PopulationPtrs& populations,
                        const PopulationDatas& population_datas)
    {
        ostringstream filename;
        filename << "genotypes_" << generation_number << ".txt"; 
        bfs::ofstream os(outdir_ / filename.str());
        if (!os)
            throw runtime_error(("[Reporter_Genotypes] Unable to open " + filename.str()).c_str());

        if (population_datas.size() != 1)
            throw runtime_error("[Reporter_Genotypes] Expecting single population.");
    
        GenotypeDataPtr genotypes = population_datas[0].genotypes->at(locus_);
        copy(genotypes->begin(), genotypes->end(), ostream_iterator<int>(os, "\n"));
    }

    virtual void update_final(size_t generation_number,
                              const PopulationPtrs& populations,
                              const PopulationDatas& population_datas)
    {}

    private:

    bfs::path outdir_;
    Locus locus_;
};


class Reporter_Fitnesses : public Reporter
{
    public:

    Reporter_Fitnesses(const string& output_directory)
    :   outdir_(output_directory)
    {}

    virtual void update(size_t generation_number,
                        const PopulationPtrs& populations,
                        const PopulationDatas& population_datas)
    {
        ostringstream filename;
        filename << "fitnesses_" << generation_number << ".txt"; 
        bfs::ofstream os(outdir_ / filename.str());
        if (!os)
            throw runtime_error(("[Reporter_Fitnesses] Unable to open " + filename.str()).c_str());

        if (population_datas.size() != 1)
            throw runtime_error("[Reporter_Fitnesses] Expecting single population.");

        const DataVector& fitnesses = *population_datas[0].fitnesses;
        copy(fitnesses.begin(), fitnesses.end(), ostream_iterator<double>(os, "\n"));
    }

    virtual void update_final(size_t generation_number,
                              const PopulationPtrs& populations,
                              const PopulationDatas& population_datas)
    {}

    private:

    bfs::path outdir_;
};


class SNPIndicator_SingleLocusHardyWeinberg : public SNPIndicator
{
    public:

    SNPIndicator_SingleLocusHardyWeinberg(Locus locus, size_t population_size, double allele_frequency)
    :   locus_(locus)
    {
        const double p = allele_frequency;
        const double q = 1-p;
        const size_t N = population_size;

        // genotypes:  (2, ... , 2, 1, ... , 1, 0, ... , 0)
        //                N*p^2       N(2pq)       N*q^2
        max_2_ = N * p * p;
        max_1_ = N * (1 - q*q);

        cout << "max_2: " << max_2_ << endl;
        cout << "max_1: " << max_1_ << endl;
    }

    virtual unsigned int operator()(unsigned int chromosome_id, const Locus& locus) const
    {
        if (locus != locus_) return 0;
        Chromosome::ID id(chromosome_id);
        if (id.individual < max_2_) return 1;
        else if (id.individual < max_1_ && id.which == 0) return 1;
        return 0; 
    }

    private:

    Locus locus_;
    size_t max_2_;
    size_t max_1_;
};


class QuantitativeTrait_SingleLocusFitness : public QuantitativeTrait
{
    public:

    QuantitativeTrait_SingleLocusFitness(int id, Locus locus, vector<double> w)
    :   QuantitativeTrait(id), locus_(locus), w_(w)
    {
        loci_.insert(locus); // to satisfy interface loci() -- TODO: revisit
    }

    virtual DataVectorPtr calculate_trait_values(GenotypeMapPtr genotypes) const
    {
        if (genotypes->size() != 1 || genotypes->count(locus_) != 1)
            throw runtime_error("[QuantitativeTrait_SingleLocusFitness] Invalid genotype map.\n");

        const GenotypeDataPtr& g = genotypes->at(locus_);

        DataVectorPtr fitnesses(new DataVector(g->size()));

        // transform {0, 1, 2} -> {w[0], w[1], w[2]}
        DataVector::iterator jt = fitnesses->begin();
        for (GenotypeData::const_iterator it=g->begin(); it!=g->end(); ++it, ++jt)
            *jt = w_[*it];

        return fitnesses;
    }

    private:

    Locus locus_;
    vector<double> w_; // relative fitnesses
};


//////////////////////////////////////////////////////////////////////////////////////////////////


SimulationController_SingleLocusSelection::Config::Config(const Parameters& parameters)
:   w(3, 1)
{
    seed = parameters.count("seed") ? atoi(parameters.at("seed").c_str()) : 0;
    output_directory = parameters.count("outdir") ? parameters.at("outdir") : "";

    population_size = parameters.count("popsize") ? atoi(parameters.at("popsize").c_str()) : 0;
    generation_count = parameters.count("gencount") ? atoi(parameters.at("gencount").c_str()) : 0;
    initial_allele_frequency = parameters.count("allelefreq") ? atof(parameters.at("allelefreq").c_str()) : 0; // atof

    cout << "allelefreq string: " << parameters.at("allelefreq") << endl;
    cout << "allelefreq: " << atof(parameters.at("allelefreq").c_str()) << endl;

    w[0] = parameters.count("w0") ? atof(parameters.at("w0").c_str()) : 1;
    w[1] = parameters.count("w1") ? atof(parameters.at("w1").c_str()) : 1;
    w[2] = parameters.count("w2") ? atof(parameters.at("w2").c_str()) : 1;
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
    cout << "initial_allele_frequency: " << config_.initial_allele_frequency << endl; // TODO: remove
    cout << "w: ";
    copy(config_.w.begin(), config_.w.end(), ostream_iterator<double>(cout, " "));
    cout << endl;

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

    Locus locus(0, 100000); // hardcoded locus
    int qtid = 0; // hardcoded QT id

    QuantitativeTraitPtr qt(new QuantitativeTrait_SingleLocusFitness(qtid, locus, config_.w));
    simulator_config_.quantitative_traits.push_back(qt);

    FitnessFunctionPtr ff(new FitnessFunction_Identity(qtid));
    simulator_config_.fitness_function = ff;

    simulator_config_.snp_indicator = SNPIndicatorPtr(new SNPIndicator_SingleLocusHardyWeinberg(
        locus, config_.population_size, config_.initial_allele_frequency));
    
    simulator_config_.reporters.push_back(ReporterPtr(new Reporter_Temp(config_.output_directory)));
    simulator_config_.reporters.push_back(ReporterPtr(new Reporter_Genotypes(config_.output_directory, locus)));
    simulator_config_.reporters.push_back(ReporterPtr(new Reporter_Fitnesses(config_.output_directory)));

    // Reporters:  
    //   full population for debugging 
    //   allele freq / homozygosity
    //   block lengths?
    //   mean fitness

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



