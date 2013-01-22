//
// recombine_data.cpp
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
#include "Chromosome.hpp"
#include "MSFormat.hpp"
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
    string populationFilename;
    string msFilename1;
    string msFilename2;
    bfs::path outputDirectory;
    size_t chromosomePairIndex;
    unsigned int positionBegin;
    unsigned int positionEnd;

    Config()
    :   chromosomePairIndex(0), positionBegin(0), positionEnd(0)
    {}
};


struct RelativePosition
{
    RelativePosition(unsigned int begin, unsigned int end)
    :   begin_(begin), end_(end)
    {
        if (begin > end)
            throw runtime_error("[RelativePosition] I am insane.");
    }

    double operator()(unsigned int position) const
    {
        if (position <= begin_) return 0;
        if (position >= end_) return 1;
        return double(position-begin_)/(end_-begin_);
    }

    private:
    unsigned int begin_;
    unsigned int end_;
};


string getSequence(const Chromosome& chromosome, const MSFormat& ms1, 
                   const MSFormat& ms2, const RelativePosition& relativePosition)
{
    string sequence;

    for (vector<DNABlock>::const_iterator it=chromosome.blocks().begin(); 
         it!=chromosome.blocks().end(); ++it)
    {
        double begin = relativePosition(it->position);
        double end = (it+1 == chromosome.blocks().end()) ? 1 : relativePosition((it+1)->position);

        Chromosome::ID id(it->id);
        const MSFormat* ms = id.population == 1 ? &ms1 : &ms2;
        
        sequence += ms->sequence(id.individual, begin, end);

        // cout << id << " " << it->position << " " << begin << " " << end << endl;
    }

    return sequence;
}


void recombineData(const Config& config)
{
    Random random;

    cout << "Reading population data.\n";
    //Population population(config.populationFilename);
    Population population;
    ifstream is(config.populationFilename.c_str());
    if (!is) throw runtime_error(("[Population] Unable to open file " + config.populationFilename).c_str());
    is >> population;
    if (population.organisms().empty())
        cerr << "[Population] Warning: no population data read from file " << config.populationFilename << endl;

    cout << "Population size: " << population.organisms().size() << endl;

    MSFormat ms1(config.msFilename1);
    MSFormat ms2(config.msFilename2);

    MSFormat result;
    result.positions = ms1.positions;

    MSFormat check;
    check.positions = ms2.positions;

    if (result != check)
        throw runtime_error("[recombineData] Segsites don't match.");

    RelativePosition relativePosition(config.positionBegin, config.positionEnd);

    for (vector<Organism>::const_iterator it=population.organisms().begin(); 
         it!=population.organisms().end(); ++it)
    {
        if (config.chromosomePairIndex >= it->chromosomePairs().size())
            throw runtime_error("[recombineData()] Bad chromosomePairIndex.");

        const ChromosomePair& cp = it->chromosomePairs()[config.chromosomePairIndex];
        
        result.sequences.push_back(getSequence(cp.first, ms1, ms2, relativePosition));
        result.sequences.push_back(getSequence(cp.second, ms1, ms2, relativePosition));
    }

    bfs::ofstream os(config.outputDirectory / "sequences.txt");
    os << result;
}


void checkExistence(const string& filename)
{
    if (!bfs::exists(filename))
    {
        ostringstream oss;
        oss << "File not found: " << filename;
        throw runtime_error(oss.str());
    }
}


Config parseCommandLine(int argc, char* argv[])
{
    if (argc != 8)
    {
        cout << "Usage: recombine_data <populationFile> <msFile1> <msFile2> <outdir> <chromosome> <begin> <end>\n";
        cout << "\n";
        cout << "  populationFile :  population file output from simrecomb\n";
        cout << "  msFile1        :  ancestral population 1 in ms format\n";
        cout << "  msFile2        :  ancestral population 2 in ms format\n";
        cout << "  chromosome     :  chromosome pair index (0-based)\n";
        cout << "  begin          :  chromosome position mapped to 0\n";
        cout << "  end            :  chromosome position mapped to 1\n";
        cout << "\n";
        cout << "Darren Kessner\n";
        cout << "John Novembre Lab, UCLA\n";
        throw runtime_error("");
    }

    Config config;
    config.populationFilename = argv[1];
    config.msFilename1 = argv[2];
    config.msFilename2 = argv[3];
    config.outputDirectory = argv[4];
    config.chromosomePairIndex = atoi(argv[5]);
    config.positionBegin = atoi(argv[6]);
    config.positionEnd = atoi(argv[7]);

    checkExistence(config.populationFilename);
    checkExistence(config.msFilename1);
    checkExistence(config.msFilename2);

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
        recombineData(config);
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


