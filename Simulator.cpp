//
// Simulator.cpp
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


#include "Simulator.hpp"
#include <iostream>
#include <iterator>
#include <sstream>
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"


using namespace std;
namespace bfs = boost::filesystem;


//
// Simulator::Config
//


ostream& operator<<(ostream& os, const Simulator::Config& config)
{
    os << "seed " << config.seed << endl << endl;

    os << "geneticMapFilenames " << config.geneticMapFilenames.size() << endl; 
    for (vector<string>::const_iterator it=config.geneticMapFilenames.begin();
         it!=config.geneticMapFilenames.end(); ++it)
        os << "geneticMapFilename " << *it << endl;
    os << endl;

    for (size_t i=0; i<config.populationConfigs.size(); i++)
    {
        os << "generation " << i << endl;
        copy(config.populationConfigs[i].begin(), config.populationConfigs[i].end(), 
             ostream_iterator<Population::Config>(os,"\n"));
        os << endl;
    }

    return os;
}


istream& operator>>(istream& is, Simulator::Config& config)
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
            config.seed = atoi(tokens[1].c_str());
        else if (tokens[0] == "geneticMapFilenames")
            continue;
        else if (tokens[0] == "geneticMapFilename" && tokens.size()==2)
            config.geneticMapFilenames.push_back(tokens[1]);
        else if (tokens[0] == "generation")
            config.populationConfigs.push_back(vector<Population::Config>());
        else if (tokens[0] == "population")
        {
            config.populationConfigs.back().push_back(Population::Config());
            istringstream temp(buffer);
            temp >> config.populationConfigs.back().back();
        }
        else
            cerr << "Ignoring invalid configuration line:\n" << buffer << endl;
    }

    return is;
}


//
// Simulator
//


Simulator::Simulator()
{}


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

void initializeRecombinationMaps(const Simulator::Config& simulationConfig, const Random& random)
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
        {
            PopulationPtr p(new Population);
            p->create_organisms(*it);
            result->push_back(p);
            //result->push_back(shared_ptr<Population>(new Population(*it)));

        }
        else
        {
            PopulationPtr p(new Population);
            p->create_organisms(*it, *current, fitnesses, random);
            result->push_back(p);
            //result->push_back(shared_ptr<Population>(new Population(*it, *current, fitnesses, random)));
        }
    }

    return result; 
}


void Simulator::simulate(const Simulator::Config& simulationConfig, const string& outputDirectoryName)
{
    bfs::path outputDirectory(outputDirectoryName);

    bfs::ofstream os_config(outputDirectory / "simrecomb_config.txt");
    os_config  << simulationConfig;
    os_config.close();

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



