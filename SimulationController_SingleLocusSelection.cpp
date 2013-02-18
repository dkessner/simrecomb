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
#include "boost/filesystem/fstream.hpp"
#include "boost/lambda/lambda.hpp"
#include <iostream>
#include <sstream>


using namespace std;
namespace bfs = boost::filesystem;
using namespace boost::lambda;


//////////////////////////////////////////////////////////////////////////////////////////////////
//
// TODO:  move this stuff


class Reporter_Population : public Reporter
{
    public:

    Reporter_Population(const string& output_directory, bool verbose)
    :   outdir_(output_directory), verbose_(verbose)
    {}

    virtual void update(size_t generation_number,
                        const PopulationPtrs& populations,
                        const PopulationDatas& population_datas)
    {
        if (!verbose_) return;

        if (populations.size() != population_datas.size())
            throw runtime_error("[Reporter_Population] Population data size mismatch.");

        const size_t population_count = populations.size();
        const char* filestem = "population";

        for (size_t population_index=0; population_index<population_count; ++population_index)
        {
            ostringstream filename;
            filename << filestem << "_" << generation_number << "_" << population_index << ".txt"; 

            bfs::ofstream os(outdir_ / filename.str());
            if (!os)
                throw runtime_error(("[Reporter_Population] Unable to open " + filename.str()).c_str());

            os << *populations[population_index];
            os.close();
        }
    }

    virtual void update_final(size_t generation_number,
                              const PopulationPtrs& populations,
                              const PopulationDatas& population_datas)
    {}

    private:

    bfs::path outdir_;
    bool verbose_;
};


class Reporter_Genotypes : public Reporter
{
    public:

    Reporter_Genotypes(const string& output_directory, Locus locus, bool verbose)
    :   outdir_(output_directory), locus_(locus), verbose_(verbose)
    {
        os_allele_freqs_.open(outdir_ / "allele_freqs.txt");
        if (!os_allele_freqs_)
            throw runtime_error("[Reporter_Genotypes] Unable to open file allele_freqs.txt");
    }

    virtual void update(size_t generation_number,
                        const PopulationPtrs& populations,
                        const PopulationDatas& population_datas)
    {
        if (populations.size() != population_datas.size())
            throw runtime_error("[Reporter_Genotypes] Population data size mismatch.");

        const size_t population_count = populations.size();

        // update allele frequenices

        for (size_t population_index=0; population_index<population_count; ++population_index)
        {
            GenotypeDataPtr genotypes = population_datas[population_index].genotypes->at(locus_);
            os_allele_freqs_ << genotypes->allele_frequency() << " ";
        }
        os_allele_freqs_ << endl;

        // if verbose, write full genotype data

        if (!verbose_) return;

        const char* filestem = "genotypes";

        for (size_t population_index=0; population_index<population_count; ++population_index)
        {
            ostringstream filename;
            filename << filestem << "_" << generation_number << "_" << population_index << ".txt"; 

            bfs::ofstream os(outdir_ / filename.str());
            if (!os)
                throw runtime_error(("[Reporter_Genotypes] Unable to open " + filename.str()).c_str());

            GenotypeDataPtr genotypes = population_datas[population_index].genotypes->at(locus_);
            copy(genotypes->begin(), genotypes->end(), ostream_iterator<int>(os, "\n"));
            os.close();
        }
    }

    virtual void update_final(size_t generation_number,
                              const PopulationPtrs& populations,
                              const PopulationDatas& population_datas)
    {}

    private:

    bfs::path outdir_;
    bfs::ofstream os_allele_freqs_;
    Locus locus_;
    bool verbose_;
};


class Reporter_Fitnesses : public Reporter
{
    public:

    Reporter_Fitnesses(const string& output_directory, bool verbose)
    :   outdir_(output_directory), verbose_(verbose)
    {
        os_mean_.open(outdir_ / "mean_fitnesses.txt");
        if (!os_mean_)
            throw runtime_error("[Reporter_Fitnesses] Unable to open file mean_fitnesses.txt");
    }

    virtual void update(size_t generation_number,
                        const PopulationPtrs& populations,
                        const PopulationDatas& population_datas)
    {
        if (populations.size() != population_datas.size())
            throw runtime_error("[Reporter_Fitnesses] Population data size mismatch.");

        const size_t population_count = populations.size();

        // update mean fitnesses

        for (size_t population_index=0; population_index<population_count; ++population_index)
        {
            const DataVector& fitnesses = *population_datas[population_index].fitnesses;
            os_mean_ << fitnesses.mean() << " ";
        }
        os_mean_ << endl;

        // if verbose, report full fitnesses

        if (!verbose_) return;

        const char* filestem = "fitnesses";

        for (size_t population_index=0; population_index<population_count; ++population_index)
        {
            ostringstream filename;
            filename << filestem << "_" << generation_number << "_" << population_index << ".txt"; 

            bfs::ofstream os(outdir_ / filename.str());
            if (!os)
                throw runtime_error(("[Reporter_Fitnesses] Unable to open " + filename.str()).c_str());

            const DataVector& fitnesses = *population_datas[population_index].fitnesses;
            copy(fitnesses.begin(), fitnesses.end(), ostream_iterator<double>(os, "\n"));
            os.close();
        }
    }

    virtual void update_final(size_t generation_number,
                              const PopulationPtrs& populations,
                              const PopulationDatas& population_datas)
    {}

    private:

    bfs::path outdir_;
    bfs::ofstream os_mean_;
    bool verbose_;
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

        // individuals: (0, 1, 2, ... ,     max_2, ... ,     max_1, ...    )
        // genotypes:   (2, 2, 2, ... , 2,  1, 1,  ... , 1,  0, 0,  ... , 0)
        //                       N*p^2            N(2pq)          N*q^2

        max_2_ = N * p * p;
        max_1_ = N * (1 - q*q);
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

    population_count = parameters.count("popcount") ? atoi(parameters.at("popcount").c_str()) : 1;
    population_size = parameters.count("popsize") ? atoi(parameters.at("popsize").c_str()) : 0;
    generation_count = parameters.count("gencount") ? atoi(parameters.at("gencount").c_str()) : 0;
    initial_allele_frequency = parameters.count("allelefreq") ? atof(parameters.at("allelefreq").c_str()) : 0; // atof

    w[0] = parameters.count("w0") ? atof(parameters.at("w0").c_str()) : 1;
    w[1] = parameters.count("w1") ? atof(parameters.at("w1").c_str()) : 1;
    w[2] = parameters.count("w2") ? atof(parameters.at("w2").c_str()) : 1;

    verbose = parameters.count("verbose"); // no good for verbose=0
}


std::ostream& operator<<(std::ostream& os, const SimulationController_SingleLocusSelection::Config& config)
{
    os << "seed = " << config.seed << endl;
    os << "outdir = " << config.output_directory << endl;
    os << "popcount = " << config.population_count << endl;
    os << "popsize = " << config.population_size << endl;
    os << "gencount = " << config.generation_count << endl;
    os << "allelefreq = " << config.initial_allele_frequency << endl;
    os << "w0 = " << config.w[0] << endl;
    os << "w1 = " << config.w[1] << endl;
    os << "w2 = " << config.w[2] << endl;
    if (config.verbose) os << "verbose" << endl;
    return os;
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
    cout << "  popcount=<population_count>              (default: 1)\n";
    cout << "  allelefreq=<initial_allele_frequency>    (default: allelefreq=0)\n";
    cout << "  w0=<relative_fitness_genotype_0>         (default: w0=1)\n";
    cout << "  w1=<relative_fitness_genotype_1>         (default: w1=1)\n";
    cout << "  w2=<relative_fitness_genotype_2>         (default: w2=1)\n";
    cout << "  verbose\n"; 
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
    bfs::ofstream os_config(bfs::path(config_.output_directory) / "config.txt");
    os_config << config_ << endl;
    os_config.close();

    cout << config_ << endl;

    // initialize simulator

    simulator_config_.seed = config_.seed;
    simulator_config_.output_directory = config_.output_directory;    

    // population configs

    Population::Configs configs_gen_0(config_.population_count);
    for (size_t i=0; i<config_.population_count; ++i)
    {
        Population::Config& popconfig = configs_gen_0[i];
        popconfig.size = config_.population_size;
        popconfig.populationID = i;
        popconfig.chromosomePairCount = 1;
    }
    simulator_config_.population_configs.push_back(configs_gen_0); // copy

    Population::Configs configs_gen_next(config_.population_count);
    for (size_t i=0; i<config_.population_count; ++i)
    {
        Population::Config& config_gen_next = configs_gen_next[i];
        config_gen_next.size = config_.population_size;
        config_gen_next.populationID = i; 
        config_gen_next.matingDistribution.push_back(1, make_pair(i,i));
    }
    for (size_t i=0; i<config_.generation_count; ++i)
        simulator_config_.population_configs.push_back(configs_gen_next); // copy

    bfs::ofstream os(bfs::path(config_.output_directory) / "popconfig.txt");
    if (!os)
        throw runtime_error("[SimulationController_SingleLocusSelection] Unable to open popconfig.txt");
    os << simulator_config_.population_configs;
    os.close();

    // QuantitativeTraits, FitnessFunction, SNPIndicator, Reporters

    Locus locus(0, 100000); // hardcoded locus
    int qtid = 0; // hardcoded QT id

    QuantitativeTraitPtr qt(new QuantitativeTrait_SingleLocusFitness(qtid, locus, config_.w));
    simulator_config_.quantitative_traits.push_back(qt);

    FitnessFunctionPtr ff(new FitnessFunction_Identity(qtid));
    simulator_config_.fitness_function = ff;

    simulator_config_.snp_indicator = SNPIndicatorPtr(new SNPIndicator_SingleLocusHardyWeinberg(
        locus, config_.population_size, config_.initial_allele_frequency));
    
    simulator_config_.reporters.push_back(ReporterPtr(new Reporter_Population(config_.output_directory, config_.verbose)));
    simulator_config_.reporters.push_back(ReporterPtr(new Reporter_Genotypes(config_.output_directory, locus, config_.verbose)));
    simulator_config_.reporters.push_back(ReporterPtr(new Reporter_Fitnesses(config_.output_directory, config_.verbose)));

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



