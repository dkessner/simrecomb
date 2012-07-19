//
// PopulationTest.cpp
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
#include "unit.hpp"
#include <iostream>
#include <cstring>


using namespace std;
using boost::shared_ptr;


ostream* os_ = 0;


void testMatingDistribution()
{
    if (os_) *os_ << "testMatingDistribution()\n";

    MatingDistribution md;
    md.push_back(.6, make_pair(0,0));
    md.push_back(.2, make_pair(1,0));
    md.push_back(.1, make_pair(2,0));

    unit_assert(md.entries().size() == 3);
    unit_assert(md.entries()[0].cumulativeWeight == .6);
    unit_assert(md.entries()[1].cumulativeWeight == .8);
    unit_assert(md.entries()[2].cumulativeWeight == .9);
    unit_assert(md.entries()[2].indexPair.first == 2);
    unit_assert(md.entries()[2].indexPair.second == 0);

    vector<int> counts(3);

    for (int i=0; i<9000; i++)
    {
        const MatingDistribution::IndexPair& p = md.random();
        counts[p.first]++;
    }

    for (size_t i=0; i<3; i++)
        if (os_) *os_ << "count " << i << ": " << counts[i] << endl;


    // test I/O

    ostringstream oss;
    oss << md;
    if (os_) *os_ << oss.str() << endl;
    MatingDistribution md2;
    md2.push_back(.42, make_pair(4,5));
    istringstream iss(oss.str());
    iss >> md2;
    unit_assert(md==md2);

    if (os_) *os_ << endl;
}


void testPopulationConfigIO()
{
    if (os_) *os_ << "testPopulationConfigIO()\n";

    Population::Config config1;
    config1.size = 4;
    config1.populationID = 1;
    config1.idOffset = 1000; // individuals==1000,1001,..., 
    config1.chromosomePairCount = 3;
    config1.matingDistribution.push_back(.42, make_pair(0,0));
    config1.matingDistribution.push_back(.66, make_pair(1,0));
    config1.matingDistribution.push_back(.23, make_pair(2,1));

    ostringstream oss;
    oss << config1;

    if (os_) *os_ << "config1:\n" << config1 << endl;

    Population::Config config2;
    unit_assert(config1 != config2);
    istringstream iss(oss.str());
    iss >> config2;

    if (os_) *os_ << "config2:\n" << config2 << endl;
    unit_assert(config1 == config2);

    if (os_) *os_ << endl;
}


void testPopulation_initial()
{
    if (os_) *os_ << "testPopulation_initial()\n";

    Population::Config config0;
    config0.size = 10;
    config0.populationID = 0;
    Population p0(config0);
    unit_assert(p0.organisms().size() == 10);
    if (os_) *os_ << "p0:\n" << p0 << endl;

    Population::Config config1;
    config1.size = 4;
    config1.populationID = 1;
    config1.idOffset = 1000; // individuals==1000,1001,..., 
    config1.chromosomePairCount = 3;
    Population p1(config1);
    unit_assert(p1.organisms().size() == 4);
    if (os_) *os_ << "p1:\n" << p1 << endl;
    if (os_) *os_ << endl;
}


void testPopulation_generated()
{
    if (os_) *os_ << "testPopulation_generated()\n";

    Population::Config config0;
    config0.size = 10;
    config0.populationID = 0;

    Population::Config config1;
    config1.size = 10;
    config1.populationID = 1;

    shared_ptr<Population> p0(new Population(config0));
    shared_ptr<Population> p1(new Population(config1));

    if (os_) *os_ << "p0:\n" << *p0 << endl;
    if (os_) *os_ << "p1:\n" << *p1 << endl;

    Populations populations;
    populations.push_back(p0);
    populations.push_back(p1);

    
    Population::Config config_nextgen;
    config_nextgen.size = 100;
    config_nextgen.matingDistribution.push_back(.5, make_pair(0,0));
    config_nextgen.matingDistribution.push_back(.3, make_pair(1,0));
    config_nextgen.matingDistribution.push_back(.2, make_pair(1,1));


    Population nextgen(config_nextgen, populations);
    if (os_) *os_ << "nextgen:\n" << nextgen << endl;

    int count00 = 0;
    int count01 = 0;
    int count11 = 0;
    
    for (vector<Organism>::const_iterator it=nextgen.organisms().begin(); it!=nextgen.organisms().end(); ++it)
    {
        Chromosome::ID id1(it->chromosomePairs()[0].first.blocks()[0].id);
        Chromosome::ID id2(it->chromosomePairs()[0].second.blocks()[0].id);

        if (id1.population==0 && id2.population==0) count00++;
        else if (id1.population==1 && id2.population==1) count11++;
        else count01++;
    }

    if (os_) *os_ << "count00: " << count00 << endl;
    if (os_) *os_ << "count01: " << count01 << endl;
    if (os_) *os_ << "count11: " << count11 << endl;
    if (os_) *os_ << endl;
}


void testPopulationIO()
{
    if (os_) *os_ << "testPopulationIO()\n";

    Population::Config config;
    config.size = 10;
    config.populationID = 0;
    Population p(config);
    unit_assert(p.organisms().size() == 10);
    if (os_) *os_ << "p:\n" << p << endl;

    ostringstream oss;
    oss << p;

    Population::Config config2;
    Population q(config2);        
    unit_assert(p != q);

    istringstream iss(oss.str());
    iss >> q;
    if (os_) *os_ << "q:\n" << q << endl;
    unit_assert(p == q);

    if (os_) *os_ << endl;
}


void test()
{
    testMatingDistribution();
    testPopulation_initial();
    testPopulation_generated();
    testPopulationConfigIO();
    testPopulationIO();
}


int main(int argc, char* argv[])
{
    try
    {
        if (argc>1 && !strcmp(argv[1],"-v")) os_ = &cout;
        test();
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


