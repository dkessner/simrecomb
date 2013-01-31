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


#include "Simulator.hpp"
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


void initializeExampleSimulatorConfig(Simulator::Config& config)
{
    const size_t chromosomePairCount_ = 3;
    const unsigned int populationSize_ = 10000;
    const double admixtureProportion_ = .8; // fraction of genes from 1st population

    for (size_t i=0; i<chromosomePairCount_; i++) 
        config.geneticMapFilenames.push_back("genetic_map_chr21_b36.txt"); // hack

    // generation 0 (ancestral populations)

    config.populationConfigs.push_back(vector<Population::Config>(3));

    Population::Config* config_pop = &config.populationConfigs[0][0];
    config_pop->size = 0;
    config_pop->populationID = 0;

    config_pop = &config.populationConfigs[0][1];
    config_pop->size = populationSize_;
    config_pop->populationID = 1;
    config_pop->chromosomePairCount = chromosomePairCount_;

    config_pop = &config.populationConfigs[0][2];
    config_pop->size = populationSize_;
    config_pop->populationID = 2;
    config_pop->chromosomePairCount = chromosomePairCount_;

    // generation 1 (initial admixture)

    config.populationConfigs.push_back(vector<Population::Config>(1));

    config_pop = &config.populationConfigs[1][0];
    config_pop->size = populationSize_;
    double p = admixtureProportion_;
    config_pop->matingDistribution.push_back(p*p, make_pair(1,1));
    config_pop->matingDistribution.push_back(2*p*(1-p), make_pair(1,2));
    config_pop->matingDistribution.push_back((1-p)*(1-p), make_pair(2,2));

    // subsequent generations - just recombination

    for (size_t generation=2; generation<8; generation++)
    {
        config.populationConfigs.push_back(vector<Population::Config>(1));
        config_pop = &config.populationConfigs[generation][0];
        config_pop->size = populationSize_;
        config_pop->matingDistribution.push_back(1, make_pair(0,0));
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
            Simulator::Config config;
            is >> config;
            is.close();

            bfs::create_directories(outputDirectory);

            Simulator simulator(config, outputDirectory);
            simulator.simulate_all();
        }
        else if (subfunction == "example")
        {
            Simulator::Config temp;
            initializeExampleSimulatorConfig(temp);
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


