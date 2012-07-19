//
// ChromosomeTest.cpp
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

#include "Chromosome.hpp"
#include "unit.hpp"
#include <iostream>
#include <map>
#include <cstring>
#include <algorithm>


using namespace std;


ostream* os_ = 0;


void test_DNABlock() 
{
    if (os_) *os_ << "test_DNABlock()\n";

    DNABlock a(1000,420);
    DNABlock b(1000,420);
    unit_assert(a==b);
    b.id = 421;
    unit_assert(a!=b);
    unit_assert(a<b);
    a.position = 1001;
    unit_assert(a>b);

    map<DNABlock,int> m;
    m[a] = 5;
    m[b] = 6;
    
    for (map<DNABlock,int>::const_iterator it=m.begin(); it!=m.end(); ++it)
        if (os_) *os_ << it->first << ": " << it->second << endl; 

    if (os_) *os_ << endl;
}


void test_DNABlock_write_read() 
{
    if (os_) *os_ << "test_DNABlock_write_read()\n";

    DNABlock a(1000,420);
    if (os_) *os_ << "a: " << a << endl;

    ostringstream oss;
    oss << a;

    istringstream iss(oss.str());
    DNABlock b(0,0);
    unit_assert(a != b);

    iss >> b;
    if (os_) *os_ << "b: " << b << endl;
    unit_assert(a == b);

    if (os_) *os_ << endl;
}


void test_id()
{
    if (os_) *os_ << "test_id()\n";

    unsigned int population = 5;
    unsigned int individual = 12345;
    unsigned int pair = 21;
    unsigned int which = 1;

    Chromosome::ID id(population, individual, pair, which);
    if (os_) *os_ << "id: " << id << endl;

    unsigned int encoded = (unsigned int)(id);
    if (os_) *os_ << "encoded: " << hex << encoded << dec << endl;

    Chromosome::ID id2(encoded);
    unit_assert(id2.population == population);
    unit_assert(id2.individual == individual);
    unit_assert(id2.pair == pair);
    unit_assert(id2.which == which);

    unit_assert(id == id2); // from automatic conversion

    if (os_) *os_ << endl;
}


void test_id_write_read()
{
    if (os_) *os_ << "test_id_write_read()\n";

    unsigned int population = 5;
    unsigned int individual = 12345;
    unsigned int pair = 21;
    unsigned int which = 1;

    Chromosome::ID id(population, individual, pair, which);
    if (os_) *os_ << "id: " << id << endl;

    ostringstream oss;
    oss << id;

    istringstream iss(oss.str());
    Chromosome::ID id2(0);
    unit_assert(id != id2);

    iss >> id2;
    unit_assert(id == id2);

    if (os_) *os_ << endl;
}


void test_construction()
{
    if (os_) *os_ << "test_construction()\n";
    Chromosome c(666);
    unit_assert(c.blocks().size() == 1);
    unit_assert(c.blocks().front().id == 666);
    if (os_) *os_ << c << endl;
    if (os_) *os_ << endl;
}


void test_equality()
{
    if (os_) *os_ << "test_equality()\n";

    DNABlocks blocks;
    for (unsigned int i=0; i<10; i++)
        blocks.push_back(DNABlock(i*1000, i));

    Chromosome a(blocks);
    Chromosome b(blocks);
    
    unit_assert(a == b);

    blocks.push_back(DNABlock(420,11));
    Chromosome c(blocks);
    unit_assert(a != c);

    if (os_) *os_ << endl;
}


void test_write_read()
{
    if (os_) *os_ << "test_write_read()\n";

    DNABlocks blocks;
    for (unsigned int i=0; i<10; i++)
        blocks.push_back(DNABlock(i*1000, i));

    Chromosome a(blocks);
    if (os_) *os_ << a << endl;

    ostringstream oss;
    oss << a;

    istringstream iss(oss.str());
    Chromosome b(0);
    unit_assert(a != b);

    iss >> b;
    if (os_) *os_ << b << endl;
    unit_assert(a == b);

    if (os_) *os_ << endl;
}


void test_extract_blocks()
{
    if (os_) *os_ << "test_extract_blocks()\n";

    unit_assert(DNABlock(1000,3) < DNABlock(2000,1)); // operator<

    DNABlocks blocks;
    for (unsigned int i=0; i<10; i++)
        blocks.push_back(DNABlock(i*1000, i));

    Chromosome c(blocks);
    if (os_) *os_ << "c: " << c << endl;

    DNABlocks extracted1;
    c.extract_blocks(666, 3666, extracted1);
    if (os_) *os_ << "extracted1: " << extracted1 << endl;
    unit_assert(extracted1.size() == 4);

    DNABlocks extracted2;
    c.extract_blocks(11000, 100000, extracted2);
    if (os_) *os_ << "extracted2: " << extracted2 << endl;
    unit_assert(extracted2.size() == 1);

    DNABlocks extracted3;
    c.extract_blocks(8500, 100000, extracted3);
    if (os_) *os_ << "extracted3: " << extracted3 << endl;
    unit_assert(extracted3.size() == 2);

    if (os_) *os_ << endl;
}


void test_recombine()
{
    if (os_) *os_ << "test_recombine()\n";
    
    Chromosome x(1);
    Chromosome y(2);

    if (os_) *os_ << "x: " << x << endl;
    if (os_) *os_ << "y: " << y << endl;

    vector<unsigned int> a_positions;
    Chromosome a(x, y, a_positions);
    if (os_) *os_ << "a: " << a << endl;
    unit_assert(a.blocks().size() == 1);
    unit_assert(a.blocks().back().id == 1);

    vector<unsigned int> b_positions;
    b_positions.push_back(0);
    Chromosome b(x, y, b_positions);
    if (os_) *os_ << "b: " << b << endl;
    unit_assert(b.blocks().size() == 1);
    unit_assert(b.blocks().back().id == 2);

    vector<unsigned int> c_positions;
    c_positions.push_back(1000);
    Chromosome c(x, y, c_positions);
    if (os_) *os_ << "c: " << c << endl;
    unit_assert(c.blocks().size() == 2);
    unit_assert(c.blocks().front().id == 1);
    unit_assert(c.blocks().back().id == 2);
    unit_assert(c.blocks().back().position == 1000);

    vector<unsigned int> d_positions;
    d_positions.push_back(0);
    d_positions.push_back(1000);
    Chromosome d(x, y, d_positions);
    if (os_) *os_ << "d: " << d << endl;
    unit_assert(d.blocks().size() == 2);
    unit_assert(d.blocks().front().id == 2);
    unit_assert(d.blocks().back().id == 1);
    unit_assert(d.blocks().back().position == 1000);

    vector<unsigned int> e_positions;
    e_positions.push_back(1000);
    e_positions.push_back(2000);
    Chromosome e(x, y, e_positions);
    if (os_) *os_ << "e: " << e << endl;
    unit_assert(e.blocks().size() == 3);
    DNABlocks::const_iterator it = e.blocks().begin();
    unit_assert(it->id == 1);
    unit_assert(it->position == 0);
    ++it;
    unit_assert(it->id == 2);
    unit_assert(it->position == 1000);
    ++it;
    unit_assert(it->id == 1);
    unit_assert(it->position == 2000);
    ++it;

    vector<unsigned int> f_positions;
    f_positions.push_back(0);
    f_positions.push_back(1000);
    f_positions.push_back(2000);
    Chromosome f(x, y, f_positions);
    if (os_) *os_ << "f: " << f << endl;
    unit_assert(f.blocks().size() == 3);
    it = f.blocks().begin();
    unit_assert(it->id == 2);
    unit_assert(it->position == 0);
    ++it;
    unit_assert(it->id == 1);
    unit_assert(it->position == 1000);
    ++it;
    unit_assert(it->id == 2);
    unit_assert(it->position == 2000);
    ++it;

    vector<unsigned int> g_positions;
    g_positions.push_back(1000);
    g_positions.push_back(2000);
    g_positions.push_back(3000);
    Chromosome g(x, y, g_positions);
    if (os_) *os_ << "g: " << g << endl;
    unit_assert(g.blocks().size() == 4);
    it = g.blocks().begin();
    unit_assert(it->id == 1);
    unit_assert(it->position == 0);
    ++it;
    unit_assert(it->id == 2);
    unit_assert(it->position == 1000);
    ++it;
    unit_assert(it->id == 1);
    unit_assert(it->position == 2000);
    ++it;
    unit_assert(it->id == 2);
    unit_assert(it->position == 3000);
    ++it;

    vector<unsigned int> h_positions;
    h_positions.push_back(0);
    h_positions.push_back(1000);
    h_positions.push_back(2000);
    h_positions.push_back(3000);
    Chromosome h(x, y, h_positions);
    if (os_) *os_ << "h: " << h << endl;
    unit_assert(h.blocks().size() == 4);
    it = h.blocks().begin();
    unit_assert(it->id == 2);
    unit_assert(it->position == 0);
    ++it;
    unit_assert(it->id == 1);
    unit_assert(it->position == 1000);
    ++it;
    unit_assert(it->id == 2);
    unit_assert(it->position == 2000);
    ++it;
    unit_assert(it->id == 1);
    unit_assert(it->position == 3000);
    ++it;

    if (os_) *os_ << endl;
}


void test_recombine_2()
{
    if (os_) *os_ << "test_recombine_2()\n";
    
    Chromosome x(1);
    Chromosome y(2);
    Chromosome z(3);
    Chromosome w(4);

    if (os_) *os_ << "x: " << x << endl;
    if (os_) *os_ << "y: " << y << endl;
    if (os_) *os_ << "z: " << z << endl;
    if (os_) *os_ << "w: " << w << endl;

    vector<unsigned int> a_positions;
    a_positions.push_back(1000);
    Chromosome a(x, y, a_positions);

    vector<unsigned int> b_positions;
    b_positions.push_back(1000);
    Chromosome b(z, w, b_positions);

    if (os_) *os_ << "a: " << a << endl;
    if (os_) *os_ << "b: " << b << endl;

    vector<unsigned int> c_positions;
    c_positions.push_back(500);
    c_positions.push_back(1500);
    c_positions.push_back(2500);
    Chromosome c(a, b, c_positions);

    if (os_) *os_ << "c: " << c << endl;

    unit_assert(c.blocks().size() == 5);
    DNABlocks::const_iterator it = c.blocks().begin();
    unit_assert(it->position == 0);
    unit_assert(it->id == 1);
    ++it;
    unit_assert(it->position == 500);
    unit_assert(it->id == 3);
    ++it;
    unit_assert(it->position == 1000);
    unit_assert(it->id == 4);
    ++it;
    unit_assert(it->position == 1500);
    unit_assert(it->id == 2);
    ++it;
    unit_assert(it->position == 2500);
    unit_assert(it->id == 4);

    vector<unsigned int> d_positions;
    d_positions.push_back(0);
    d_positions.push_back(500);
    d_positions.push_back(1500);
    d_positions.push_back(2500);
    Chromosome d(a, b, d_positions);

    if (os_) *os_ << "d: " << d << endl;

    unit_assert(d.blocks().size() == 5);
    it = d.blocks().begin();
    unit_assert(it->position == 0);
    unit_assert(it->id == 3);
    ++it;
    unit_assert(it->position == 500);
    unit_assert(it->id == 1);
    ++it;
    unit_assert(it->position == 1000);
    unit_assert(it->id == 2);
    ++it;
    unit_assert(it->position == 1500);
    unit_assert(it->id == 4);
    ++it;
    unit_assert(it->position == 2500);
    unit_assert(it->id == 2);

    if (os_) *os_ << endl;
}


void test_recombine_3()
{
    if (os_) *os_ << "test_recombine_3()\n";
    
    DNABlocks blocks_a;
    blocks_a.push_back(DNABlock(0, 1));
    blocks_a.push_back(DNABlock(1000, 2));
    blocks_a.push_back(DNABlock(2000, 3));

    DNABlocks blocks_b;
    blocks_b.push_back(DNABlock(0, 4));
    blocks_b.push_back(DNABlock(1000, 5));
    blocks_b.push_back(DNABlock(2000, 6));

    Chromosome ch_a(blocks_a);
    Chromosome ch_b(blocks_b);

    if (os_) *os_ << "ch_a: " << ch_a << endl;
    if (os_) *os_ << "ch_b: " << ch_b << endl;

    vector<unsigned int> c_positions;
    c_positions.push_back(500);
    c_positions.push_back(3000);
    c_positions.push_back(4000);
    Chromosome ch_c(ch_a, ch_b, c_positions);

    if (os_) *os_ << "ch_c: " << ch_c << endl;
            
    unit_assert(ch_c.blocks().size() == 6);
    DNABlocks::const_iterator it = ch_c.blocks().begin();
    unit_assert(it->position == 0);
    unit_assert(it->id == 1);
    ++it;
    unit_assert(it->position == 500);
    unit_assert(it->id == 4);
    ++it;
    unit_assert(it->position == 1000);
    unit_assert(it->id == 5);
    ++it;
    unit_assert(it->position == 2000);
    unit_assert(it->id == 6);
    ++it;
    unit_assert(it->position == 3000);
    unit_assert(it->id == 3);
    ++it;
    unit_assert(it->position == 4000);
    unit_assert(it->id == 6);
    
    vector<unsigned int> d_positions;
    d_positions.push_back(0);
    d_positions.push_back(500);
    d_positions.push_back(3000);
    d_positions.push_back(4000);
    Chromosome ch_d(ch_a, ch_b, d_positions);

    if (os_) *os_ << "ch_d: " << ch_d << endl;
            
    unit_assert(ch_d.blocks().size() == 6);
    it = ch_d.blocks().begin();
    unit_assert(it->position == 0);
    unit_assert(it->id == 4);
    ++it;
    unit_assert(it->position == 500);
    unit_assert(it->id == 1);
    ++it;
    unit_assert(it->position == 1000);
    unit_assert(it->id == 2);
    ++it;
    unit_assert(it->position == 2000);
    unit_assert(it->id == 3);
    ++it;
    unit_assert(it->position == 3000);
    unit_assert(it->id == 6);
    ++it;
    unit_assert(it->position == 4000);
    unit_assert(it->id == 3);

    if (os_) *os_ << endl;
}


void test_recombine_4()
{
    /*
    problem with client passing unsorted positions for recombination;
    chose not to sort internally to allow positions to be passed by const&;
    if this bites again, we should sort internally

    baby:
    + { (0,6001) (25430475,6101) }
    - { (0,101) (26347124,1) }
    + { (0,6002) }
    - { (0,102) (24256374,2) }
    + { (0,6003) }
    - { (0,103) (28440008,3) }

    baby gamete
    { (0,101) (26347124,1) (31872971,6101) (14649488,101) (26347124,1) }  <-- yikes!
    { (0,102) (24256374,2) }
    { (0,103) (28440008,3) (45267208,6003) }
    */

    if (os_) *os_ << "test_recombine_4()\n";

    DNABlocks blocks_a;
    blocks_a.push_back(DNABlock(0, 6001));
    blocks_a.push_back(DNABlock(25430475,6101));

    DNABlocks blocks_b;
    blocks_b.push_back(DNABlock(0, 101));
    blocks_b.push_back(DNABlock(26347124,1));

    Chromosome ch_a(blocks_a);
    Chromosome ch_b(blocks_b);

    if (os_) *os_ << "ch_a: " << ch_a << endl;
    if (os_) *os_ << "ch_b: " << ch_b << endl;

    vector<unsigned int> c_positions;
    c_positions.push_back(0);
    c_positions.push_back(31872971);
    c_positions.push_back(14649488);
    sort(c_positions.begin(), c_positions.end()); // sort
    Chromosome ch_c(ch_a, ch_b, c_positions);

    if (os_) *os_ << "ch_c: " << ch_c << endl;

    unit_assert(ch_c.blocks().size() == 4);
    DNABlocks::const_iterator it = ch_c.blocks().begin();
    unit_assert(it->position == 0);
    unit_assert(it->id == 101);
    ++it;
    unit_assert(it->position == 14649488);
    unit_assert(it->id == 6001);
    ++it;
    unit_assert(it->position == 25430475);
    unit_assert(it->id == 6101);
    ++it;
    unit_assert(it->position == 31872971);
    unit_assert(it->id == 1);

    if (os_) *os_ << endl;
}


void test()
{
    test_DNABlock();
    test_DNABlock_write_read();
    test_id();
    test_id_write_read();
    test_construction();
    test_equality();
    test_write_read();
    test_extract_blocks();
    test_recombine();
    test_recombine_2();
    test_recombine_3();
    test_recombine_4();
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


