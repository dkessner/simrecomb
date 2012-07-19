//
// OrganismTest.cpp
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

#include "Organism.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;
using boost::shared_ptr;


ostream* os_ = 0;


void test_construction()
{
    if (os_) *os_ << "test_construction()\n";
    Chromosome::ID id(4, 666, 0, 0);
    Organism a(id, 3); // chromosomeCount == 3
    if (os_) *os_ << "a:\n" << a << endl;

    unit_assert(a.chromosomePairs().size() == 3);
    for (size_t i=0; i<a.chromosomePairs().size(); i++)
    {
        const Chromosome& first = a.chromosomePairs()[i].first;
        const Chromosome& second = a.chromosomePairs()[i].second;

        unit_assert(first.blocks().size() == 1);
        Chromosome::ID id1(first.blocks().front().id);
        unit_assert(id1.population == 4);
        unit_assert(id1.individual == 666);
        unit_assert(id1.pair == i);
        unit_assert(id1.which == 0);

        unit_assert(second.blocks().size() == 1);
        Chromosome::ID id2(second.blocks().front().id);
        unit_assert(id2.population == 4);
        unit_assert(id2.individual == 666);
        unit_assert(id2.pair == i);
        unit_assert(id2.which == 1);
    }
}


void test_equality()
{
    if (os_) *os_ << "test_equality()\n";

    Chromosome::ID id(4, 666, 0, 0);
    Organism a(id, 3); // chromosomeCount == 3

    Organism b(id, 4);
    unit_assert(a != b);

    b = Organism(id, 3);
    unit_assert(a == b);

    Chromosome::ID id2(4, 667, 0, 0);
    b = Organism(id2, 3);
    unit_assert(a != b);

    if (os_) *os_ << endl;
}


void test_write_read()
{
    if (os_) *os_ << "test_write_read()\n";

    Chromosome::ID id(4, 666, 0, 0);
    Organism a(id, 3); // chromosomeCount == 3
    if (os_) *os_ << "a:\n" << a;

    ostringstream oss;
    oss << a;

    Organism b(0);
    unit_assert(a != b);

    istringstream iss(oss.str());
    iss >> b;
    if (os_) *os_ << "b:\n" << b;
    unit_assert(a == b);

    if (os_) *os_ << endl;
}


void test_construction_gamete()
{
    if (os_) *os_ << "test_construction_gamete()\n";

    Organism::Gamete g1;
    g1.push_back(Chromosome(1));
    g1.push_back(Chromosome(2));
    g1.push_back(Chromosome(3));

    Organism::Gamete g2;
    g2.push_back(Chromosome(101));
    g2.push_back(Chromosome(102));
    g2.push_back(Chromosome(103));

    Organism baby(g1, g2);
    if (os_) *os_ << "baby:\n" << baby << endl;
    unit_assert(baby.chromosomePairs().size() == 3);
    for (size_t i=0; i<3; i++)
    {
        unit_assert(baby.chromosomePairs()[i].first.blocks().size() == 1 &&
                    baby.chromosomePairs()[i].first.blocks().front().id == 1+i);
        unit_assert(baby.chromosomePairs()[i].second.blocks().size() == 1 &&
                    baby.chromosomePairs()[i].second.blocks().front().id == 101+i);
    }
}


void test_RecombinationPositionGenerator_Trivial()
{
    if (os_) *os_ << "test_RecombinationPositionGenerator_Trivial()\n";
    RecombinationPositionGenerator_Trivial trivial;
    RecombinationPositionGenerator& r = trivial; // default param in base interface only
    for (int i=0; i<10; i++)
        if (os_) *os_ << r.get_positions().size() << endl;
    if (os_) *os_ << endl;
}


void print_gamete(ostream& os, const Organism::Gamete& g, const string& name = "gamete")
{
    os << name << endl;
    copy(g.begin(), g.end(), ostream_iterator<Chromosome>(os, "\n"));
    os << endl;
}


void test_create_gamete()
{
    if (os_) *os_ << "test_create_gamete()\n";

    Organism::Gamete g1;
    g1.push_back(Chromosome(1));
    g1.push_back(Chromosome(2));
    g1.push_back(Chromosome(3));

    Organism::Gamete g2;
    g2.push_back(Chromosome(101));
    g2.push_back(Chromosome(102));
    g2.push_back(Chromosome(103));

    Organism o(g1, g2);
    
    if (os_) 
    {
        *os_ << "organism:\n" << o << endl;
        print_gamete(*os_, o.create_gamete());
        print_gamete(*os_, o.create_gamete());
        print_gamete(*os_, o.create_gamete());
    }
}


void test_recombination_map()
{
    if (os_) *os_ << "test_recombination_map()\n";

    vector<string> filenames;
    for (int i=0; i<3; i++) 
        filenames.push_back("genetic_map_chr21_b36.txt");
    
    Organism::recombinationPositionGenerator_ =
        shared_ptr<RecombinationPositionGenerator>(
            new RecombinationPositionGenerator_RecombinationMap(filenames));

    Organism::Gamete m1;
    m1.push_back(Chromosome(1));
    m1.push_back(Chromosome(2));
    m1.push_back(Chromosome(3));

    Organism::Gamete m2;
    m2.push_back(Chromosome(101));
    m2.push_back(Chromosome(102));
    m2.push_back(Chromosome(103));

    //Organism mom(m1, m2);
    Organism mom(Chromosome::ID(0,0,0,0), 3);
    if (os_) *os_ << "mom:\n" << mom << endl;

    Organism::Gamete d1;
    d1.push_back(Chromosome(6001));
    d1.push_back(Chromosome(6002));
    d1.push_back(Chromosome(6003));

    Organism::Gamete d2;
    d2.push_back(Chromosome(6101));
    d2.push_back(Chromosome(6102));
    d2.push_back(Chromosome(6103));

    //Organism dad(d1, d2);
    Organism dad(Chromosome::ID(0,1,0,0), 3);
    if (os_) *os_ << "dad:\n" << dad << endl;
 
    Organism::Gamete egg = mom.create_gamete();
    if (os_) print_gamete(*os_, egg, "egg");
    Organism::Gamete sperm = dad.create_gamete();
    if (os_) print_gamete(*os_, sperm, "sperm");

    Organism baby(sperm, egg);
    if (os_) *os_ << "baby:\n" << baby << endl;

    if (os_)
    {
        print_gamete(*os_, baby.create_gamete(), "baby gamete");
        print_gamete(*os_, baby.create_gamete(), "baby gamete");
        print_gamete(*os_, baby.create_gamete(), "baby gamete");
        print_gamete(*os_, baby.create_gamete(), "baby gamete");
    }
}


void test()
{
    test_construction();
    test_equality();
    test_write_read();
    test_construction_gamete();
    test_RecombinationPositionGenerator_Trivial();
    test_create_gamete();
    test_recombination_map();
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


