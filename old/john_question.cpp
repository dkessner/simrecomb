//
// john_question.cpp
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
};


unsigned int begin_ = 10000000;
unsigned int end_ = 30000000;
unsigned int step1_ = 500000;
unsigned int step2_ = 100000;
unsigned int max_distance_ = 50000000;


unsigned int popid(const Organism& organism, unsigned int position)
{
    const Chromosome& ch = organism.chromosomePairs().front().first;

    vector<DNABlock>::const_iterator it = lower_bound(ch.blocks().begin(), ch.blocks().end(), 
        DNABlock(position));

    if (it != ch.blocks().begin()) it = it-1;

    Chromosome::ID id(it->id);

    //cout << position << " " << *it << " " << id.population << endl;

    return id.population;
}


double calculate_D(const Population& p, unsigned int pos1, unsigned int pos2)
{
    unsigned int count1=0, count2=0, count=0, total=0;

    for (vector<Organism>::const_iterator it=p.organisms().begin(); it!=p.organisms().end(); ++it)
    {
        unsigned int popid1 = popid(*it, pos1);
        unsigned int popid2 = popid(*it, pos2);

        if (popid1 == 1) count1++;
        if (popid2 == 1) count2++;
        if (popid1==1 && popid2==1) count++;
        total++;
    }

    //cout << count1 << " " << count2 << " " << count << " " << total << endl;

    return double(count)/total - double(count1)*count2/total/total;
}


void go(const Config& config)
{
    Population p(config.populationFilename);

    for (size_t pos1=begin_; pos1<end_; pos1 += step1_)
    for (size_t pos2=pos1+step2_; pos2<pos1+max_distance_; pos2 += step2_)
    {
        unsigned int d = pos2 - pos1;
        double D = calculate_D(p, pos1, pos2);
        cout << d << " " << D << endl;
    }
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
    if (argc != 2)
    {
        cout << "Usage: john_question <populationFile>\n";
        cout << "\n";
        cout << "Darren Kessner\n";
        cout << "John Novembre Lab, UCLA\n";
        throw runtime_error("");
    }

    Config config;
    config.populationFilename = argv[1];

    checkExistence(config.populationFilename);

    return config;
}


int main(int argc, char* argv[])
{
    try
    {
        Config config = parseCommandLine(argc, argv); 
        go(config);
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


