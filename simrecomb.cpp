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
    unsigned int seed;
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


void initializeExampleSimulationConfig(SimulationConfig& simulationConfig)
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


void initializeRecombinationMaps(const SimulationConfig& simulationConfig, const Random& random)
{
    Organism::recombinationPositionGenerator_ =
        shared_ptr<RecombinationPositionGenerator>(
            new RecombinationPositionGenerator_RecombinationMap(simulationConfig.geneticMapFilenames, random));
}


shared_ptr<Populations> createPopulations(shared_ptr<Populations> current, 
                                          const vector<Population::Config>& populationConfigs,
                                          const DataVectorPtrs& fitnesses,
                                          const Random& random)
{
    shared_ptr<Populations> result(new Populations);

    for (vector<Population::Config>::const_iterator it=populationConfigs.begin(); it!=populationConfigs.end(); ++it)
    {
        if (!current.get())
            result->push_back(shared_ptr<Population>(new Population(*it)));
        else
            result->push_back(shared_ptr<Population>(new Population(*it, *current, fitnesses, random)));
    }

    return result; 
}


void simulate(const SimulationConfig& simulationConfig, bfs::path outputDirectory)
{
    bfs::ofstream osSimulationConfig(outputDirectory / "simrecomb_config.txt");
    osSimulationConfig << simulationConfig;
    osSimulationConfig.close();

    Random random(simulationConfig.seed);

    cout << "Initializing recombination maps.\n";
    initializeRecombinationMaps(simulationConfig, random);

    bfs::ofstream osLog(outputDirectory / "log.txt");

    shared_ptr<Populations> current;
    const size_t generation_count = simulationConfig.populationConfigs.size();
    for (size_t generation=0; generation<generation_count; generation++)
    {
        cout << "Generation " << generation << endl;

        DataVectorPtrs dummy_fitnesses(current.get() ? current->size() : 0);

        shared_ptr<Populations> next = 
            createPopulations(current, simulationConfig.populationConfigs[generation], dummy_fitnesses, random);

        current = next;

        ostringstream label;
        label << "gen" << generation;
        //printPopulations(osLog, *current, label.str());
        osLog << generation << endl;
    }

    osLog.close();

    // output the last generation
    
    for (size_t i=0; i<current->size(); ++i)
    {
        ostringstream filename;
        filename << "pop" << i << ".txt"; 
        bfs::ofstream os_pop(outputDirectory / filename.str());
        os_pop << *(*current)[i];
    }
}


int main(int argc, char* argv[])
{
    try
    {
        ostringstream usage;
        usage << "Usage: simrecomb <subfunction> [args]\n";
        usage << endl;
        usage << "Print example config file to stdout:\n";
        usage << "       simrecomb example\n";
        usage << endl;
        usage << "Run simulation:\n";
        usage << "       simrecomb sim <config_filename> <outputdir>\n";
        usage << endl;
        usage << "Darren Kessner\n";
        usage << "John Novembre Lab, UCLA\n";

        string subfunction = argc>1 ? argv[1] : "";

        if (subfunction == "sim")
        {
            if (argc < 4) throw runtime_error(usage.str().c_str());
                
            string configFilename = argv[2];
            string outputDirectory = argv[3];

            if (!bfs::exists(configFilename))
                throw runtime_error(("[simrecomb] Config file not found: " + configFilename).c_str());

            if (bfs::exists(outputDirectory))
                throw runtime_error(("[simrecomb] Directory exists: " + outputDirectory).c_str());

            cout << "Reading configuration file " << configFilename << endl;
            bfs::ifstream is(configFilename);
            SimulationConfig simulationConfig;
            is >> simulationConfig;
            is.close();

            bfs::create_directories(outputDirectory);

            simulate(simulationConfig, outputDirectory);
        }
        else if (subfunction == "example")
        {
            SimulationConfig temp;
            initializeExampleSimulationConfig(temp);
            cout << temp;
        }
        else if (subfunction == "txt2pop")
        {
            if (argc < 4) throw runtime_error(usage.str().c_str());
            string filename_in = argv[2];
            string filename_out = argv[3];

            cout << "reading " << filename_in << endl << flush;
            ifstream is(filename_in.c_str());
            if (!is) throw runtime_error(("[simrecomb] Unable to open file " + filename_in).c_str());
            Population p;
            is >> p;

            cout << "writing " << filename_out << endl << flush;
            ofstream os(filename_out.c_str(), ios::binary);
            p.write(os);
            os.close();
        }
        else if (subfunction == "pop2txt")
        {
            if (argc < 4) throw runtime_error(usage.str().c_str());
            string filename_in = argv[2];
            string filename_out = argv[3];

            cout << "reading " << filename_in << endl << flush;
            ifstream is(filename_in.c_str(), ios::binary);
            if (!is) throw runtime_error(("[simrecomb] Unable to open file " + filename_in).c_str());
            Population p;
            p.read(is);

            cout << "writing " << filename_out << endl << flush;
            ofstream os(filename_out.c_str());
            os << p;
            os.close();
        }
        else
        {
            throw runtime_error(usage.str().c_str());
        }

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


