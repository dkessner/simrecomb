//
// subsample_population.cpp
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
#include "Random.hpp"
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


struct Config
{
    bfs::path filename;
    size_t sampleSize;
    size_t replicateCount;
    bfs::path outputDirectory;

    Config()
    :   sampleSize(0), replicateCount(0)
    {}
};


void subsamplePopulation(const Config& config)
{
    Random random;

    Population::Config dummy;
    Population p(dummy);
    bfs::ifstream is(config.filename);

    cout << "Reading population data.\n";
    is >> p;
    is.close();
    if (p.organisms().empty())
        throw runtime_error("Error reading population data.");

    for (size_t i=0; i<config.replicateCount; i++)
    {
        shared_ptr<Population> subsample = p.randomSubsample(config.sampleSize, random);

        ostringstream filename;
        filename << "subsample_" << config.sampleSize << "_" << i << ".txt";
        bfs::ofstream os(config.outputDirectory / filename.str());

        os << *subsample;
        os.close();
    }
}


Config parseCommandLine(int argc, char* argv[])
{
    if (argc != 5)
    {
        cout << "Usage: subsample_population <filename> <sample_size> <replicate_count> <outputdir>\n";
        cout << "\n";
        cout << "Darren Kessner\n";
        cout << "John Novembre Lab, UCLA\n";
        throw runtime_error("");
    }

    Config config;
    config.filename = argv[1];
    config.sampleSize = atoi(argv[2]);
    config.replicateCount = atoi(argv[3]);
    config.outputDirectory = argv[4];

    if (!bfs::exists(config.filename))
    {
        ostringstream oss;
        oss << "File not found: " << config.filename;
        throw runtime_error(oss.str());
    }

    if (bfs::exists(config.outputDirectory))
    {
        ostringstream oss;
        oss << "File/directory already exists: " << config.outputDirectory;
        throw runtime_error(oss.str());
    }

    bfs::create_directories(config.outputDirectory);

    return config;
}


int main(int argc, char* argv[])
{
    try
    {
        Config config = parseCommandLine(argc, argv); 
        subsamplePopulation(config);
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


