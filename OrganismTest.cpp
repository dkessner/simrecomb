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


ostream* os_ = 0;
//ostream* os_ = &cout;


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


void demo_RecombinationPositionGenerator_Trivial()
{
    if (os_) *os_ << "demo_RecombinationPositionGenerator_Trivial()\n";
    Random random;
    RecombinationPositionGenerator_Trivial trivial(random);
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


void demo_create_gamete()
{
    if (os_) *os_ << "demo_create_gamete()\n";

    Random random;

    Organism::recombinationPositionGenerator_ =
        shared_ptr<RecombinationPositionGenerator>(
            new RecombinationPositionGenerator_Trivial(random));


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


void demo_recombination_map()
{
    if (os_) *os_ << "demo_recombination_map()\n";

    vector<string> filenames;
    for (int i=0; i<3; i++) 
        filenames.push_back("genetic_map_chr21_b36.txt");

    Random random;
    
    Organism::recombinationPositionGenerator_ =
        shared_ptr<RecombinationPositionGenerator>(
            new RecombinationPositionGenerator_RecombinationMap(filenames, random));

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


class RecombinationPositionGenerator_Testing : public RecombinationPositionGenerator
{
    public:

    RecombinationPositionGenerator_Testing()
    :   count_(1)
    {}

    virtual vector<unsigned int> get_positions(size_t index) const
    {
        vector<unsigned int> result;
        if (count_ % 2) result.push_back(0); // start with 2nd chromosome when count_ is odd
        result.push_back(count_ * 10000);
        ++count_;
        return result;
    }

    private:
    mutable size_t count_;
};


void test_construction_parents()
{
    if (os_) *os_ << "test_construction_parents()\n";

    Organism mom(Chromosome::ID(0,7,0,0), 3);
    Organism dad(Chromosome::ID(0,8,0,0), 3);

    if (os_) *os_ << "mom:\n" << mom << "dad:\n" << dad;

    Organism::recombinationPositionGenerator_ = 
        shared_ptr<RecombinationPositionGenerator>(new RecombinationPositionGenerator_Testing);

    // make baby with Organism(Organism&, Organism&) constructor

    Organism baby(mom, dad);

    if (os_) *os_ << "baby:\n" << baby << endl;

    for (size_t i=0; i<6; ++i)
    {
        const ChromosomePair& p = baby.chromosomePairs()[i/2];
        const Chromosome& c = (i%2==0) ? p.first : p.second;
        Chromosome::ID id0(c.blocks()[0].id);
        Chromosome::ID id1(c.blocks()[1].id);

        unit_assert(c.blocks().size() == 2);
        unit_assert(c.blocks()[0].position == 0);
        unit_assert(c.blocks()[1].position == (i+1)*10000);
        if (i%2 == 0) 
        {
            unit_assert(id0.individual == 7 && id1.individual == 7);
            unit_assert(id0.which == 1 && id1.which == 0);
        }
        else
        {
            unit_assert(id0.individual == 8 && id1.individual == 8);
            unit_assert(id0.which == 0 && id1.which == 1);
        }
        unit_assert(id0.pair == i/2 && id1.pair == i/2);
    }
}


void test_write_read_binary()
{
    if (os_) *os_ << "test_write_read_binary()\n";

    Chromosome::ID id(4, 666, 0, 0);
    Organism a(id, 3); // chromosomeCount == 3
    if (os_) *os_ << "a:\n" << a;

    ostringstream oss;
    a.write(oss);

    Organism b;
    unit_assert(a != b);

    istringstream iss(oss.str());
    b.read(iss);
    if (os_) *os_ << "b:\n" << b;
    unit_assert(a == b);

    // harder test

    Organism::recombinationPositionGenerator_ = 
        shared_ptr<RecombinationPositionGenerator>(new RecombinationPositionGenerator_Testing);

    Organism a2(Chromosome::ID(5, 666, 0, 0), 3); // chromosomeCount == 3
    Organism child(a, a2);

    if (os_) *os_ << "a2:\n" << a2;
    if (os_) *os_ << "child:\n" << child;
    
    ostringstream oss2;
    child.write(oss2);

    istringstream iss2(oss2.str());
    Organism child_test;
    unit_assert(child != child_test);

    child_test.read(iss2);
    if (os_) *os_ << "child_test:\n" << child_test;
    unit_assert(child == child_test);
}


void test()
{
    test_construction();
    test_equality();
    test_write_read();
    test_construction_gamete();
    demo_RecombinationPositionGenerator_Trivial();
    demo_create_gamete();
    demo_recombination_map();
    test_construction_parents();
    test_write_read_binary();
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


