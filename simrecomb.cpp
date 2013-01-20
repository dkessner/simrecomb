//
// simrecomb.cpp
//
// Copyright 2012 Darren Kessner
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

#include "Population.hpp"
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <map>
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"


using namespace std;
using boost::shared_ptr;
namespace bfs = boost::filesystem;


struct SimulationConfig
{
    double seed;
    vector<string> geneticMapFilenames;
    vector< vector<Population::Config> > populationConfigs; // each generation is a vector<Population::Config>

    SimulationConfig() : seed(0) {}
};


ostream& operator<<(ostream& os, const SimulationConfig& simulationConfig)
{
    os << "seed " << simulationConfig.seed << endl << endl;

    os << "geneticMapFilenames " << simulationConfig.geneticMapFilenames.size() << endl; 
    for (vector<string>::const_iterator it=simulationConfig.geneticMapFilenames.begin();
         it!=simulationConfig.geneticMapFilenames.end(); ++it)
        os << "geneticMapFilename " << *it << endl;
    os << endl;

    for (size_t i=0; i<simulationConfig.populationConfigs.size(); i++)
    {
        os << "generation " << i << endl;
        copy(simulationConfig.populationConfigs[i].begin(), simulationConfig.populationConfigs[i].end(), 
             ostream_iterator<Population::Config>(os,"\n"));
        os << endl;
    }

    return os;
}


istream& operator>>(istream& is, SimulationConfig& simulationConfig)
{
    while (is)
    {
        // parse line by line

        string buffer;
        getline(is, buffer);
        if (!is) return is;

        vector<string> tokens;
        istringstream iss(buffer);
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

        // switch on first token

        if (tokens.empty())
            continue;
        else if (tokens[0] == "seed" && tokens.size() == 2)
            simulationConfig.seed = atoi(tokens[1].c_str());
        else if (tokens[0] == "geneticMapFilenames")
            continue;
        else if (tokens[0] == "geneticMapFilename" && tokens.size()==2)
            simulationConfig.geneticMapFilenames.push_back(tokens[1]);
        else if (tokens[0] == "generation")
            simulationConfig.populationConfigs.push_back(vector<Population::Config>());
        else if (tokens[0] == "population")
        {
            simulationConfig.populationConfigs.back().push_back(Population::Config());
            istringstream temp(buffer);
            temp >> simulationConfig.populationConfigs.back().back();
        }
        else
            cerr << "Ignoring invalid configuration line:\n" << buffer << endl;
    }

    return is;
}


struct Config
{
    bfs::path configFilename;
    bfs::path outputDirectory;
    SimulationConfig simulationConfig;
};


void printPopulations(ostream& os, const Populations& populations, const string& label)
{
    os << "<" << label << ">\n";
    for (size_t i=0; i<populations.size(); i++)
    {
        os << "<population" << i << ">\n" 
             << *populations[i]
             << "</population" << i << ">\n";
    }
    os << "</" << label << ">\n";
}


void initializeDefaultSimulationConfig(SimulationConfig& simulationConfig)
{
    const size_t chromosomePairCount_ = 3;
    const unsigned int populationSize_ = 10000;
    const double admixtureProportion_ = .8; // fraction of genes from 1st population

    for (size_t i=0; i<chromosomePairCount_; i++) 
        simulationConfig.geneticMapFilenames.push_back("genetic_map_chr21_b36.txt"); // hack

    // generation 0 (ancestral populations)

    simulationConfig.populationConfigs.push_back(vector<Population::Config>(3));

    Population::Config* config = &simulationConfig.populationConfigs[0][0];
    config->size = 0;
    config->populationID = 0;

    config = &simulationConfig.populationConfigs[0][1];
    config->size = populationSize_;
    config->populationID = 1;
    config->chromosomePairCount = chromosomePairCount_;

    config = &simulationConfig.populationConfigs[0][2];
    config->size = populationSize_;
    config->populationID = 2;
    config->chromosomePairCount = chromosomePairCount_;

    // generation 1 (initial admixture)

    simulationConfig.populationConfigs.push_back(vector<Population::Config>(1));

    config = &simulationConfig.populationConfigs[1][0];
    config->size = populationSize_;
    double p = admixtureProportion_;
    config->matingDistribution.push_back(p*p, make_pair(1,1));
    config->matingDistribution.push_back(2*p*(1-p), make_pair(1,2));
    config->matingDistribution.push_back((1-p)*(1-p), make_pair(2,2));

    // subsequent generations - just recombination

    for (size_t generation=2; generation<8; generation++)
    {
        simulationConfig.populationConfigs.push_back(vector<Population::Config>(1));
        config = &simulationConfig.populationConfigs[generation][0];
        config->size = populationSize_;
        config->matingDistribution.push_back(1, make_pair(0,0));
    }
}


void initializeRecombinationMaps(const Config& config, const Random& random)
{
    Organism::recombinationPositionGenerator_ =
        shared_ptr<RecombinationPositionGenerator>(
            new RecombinationPositionGenerator_RecombinationMap(config.simulationConfig.geneticMapFilenames, random));
}


shared_ptr<Populations> createPopulations(shared_ptr<Populations> current, 
                                          const vector<Population::Config>& populationConfigs,
                                          const Random& random)
{
    shared_ptr<Populations> result(new Populations);

    for (vector<Population::Config>::const_iterator it=populationConfigs.begin(); it!=populationConfigs.end(); ++it)
    {
        if (!current.get())
            result->push_back(shared_ptr<Population>(new Population(*it, random)));
        else
            result->push_back(shared_ptr<Population>(new Population(*it, *current, random)));
    }

    return result; 
}


void simulate(const Config& config)
{
    bfs::ofstream osSimulationConfig(config.outputDirectory / "simrecomb_config.txt");
    osSimulationConfig << config.simulationConfig;
    osSimulationConfig.close();

    Random random(config.simulationConfig.seed);

    cout << "Initializing recombination maps.\n";
    initializeRecombinationMaps(config, random);

    bfs::ofstream osLog(config.outputDirectory / "log.txt");

    shared_ptr<Populations> current;
    for (size_t generation=0; generation<config.simulationConfig.populationConfigs.size(); generation++)
    {
        cout << "Generation " << generation << endl;
        shared_ptr<Populations> next = 
            createPopulations(current, config.simulationConfig.populationConfigs[generation], random);
        current = next;
        ostringstream label;
        label << "gen" << generation;
        printPopulations(osLog, *current, label.str());
    }

    osLog.close();

    // output the last generation
    
    for (size_t i=0; i<current->size(); ++i)
    {
        ostringstream filename;
        filename << "pop" << i << ".txt"; 
        bfs::ofstream os_pop(config.outputDirectory / filename.str());
        os_pop << *(*current)[i];
    }
}


Config parseCommandLine(int argc, char* argv[])
{
    if (argc==2 && !strcmp(argv[1],"default"))
    {
        SimulationConfig temp;
        initializeDefaultSimulationConfig(temp);
        cout << temp;
        exit(0);
    }

    if (argc != 3)
    {
        cout << "Usage: simrecomb <config_filename> <outputdir>\n";
        cout << "   or: simrecomb default (prints default config to stdout)\n";
        cout << "\n";
        cout << "Darren Kessner\n";
        cout << "John Novembre Lab, UCLA\n";
        throw runtime_error("");
    }

    Config config;
    config.configFilename = argv[1];
    config.outputDirectory = argv[2];

    if (!bfs::exists(config.configFilename))
    {
        ostringstream oss;
        oss << "Config file not found: " << config.configFilename;
        throw runtime_error(oss.str());
    }

    if (bfs::exists(config.outputDirectory))
    {
        ostringstream oss;
        oss << "File/directory already exists: " << config.outputDirectory;
        throw runtime_error(oss.str());
    }

    cout << "Reading configuration file " << config.configFilename << endl;
    bfs::ifstream is(config.configFilename);
    is >> config.simulationConfig;
    is.close();

    bfs::create_directories(config.outputDirectory);

    return config;
}


int main(int argc, char* argv[])
{
    try
    {
        Config config = parseCommandLine(argc, argv); 
        simulate(config);
        return 0;
    }
    catch(exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }
    catch(...)
    {
        cerr << "Caught unknown exception.\n";
        return 1;
    }
}


