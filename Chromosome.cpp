//
// Chromosome.cpp
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
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <sstream>


using namespace std;


//
// DNABlock
//


ostream& operator<<(ostream& os, const DNABlock& x)
{
    os << "(" << x.position << "," << Chromosome::ID(x.id) << ")";
    return os;
}


istream& operator>>(istream& is, DNABlock& x)
{
    Chromosome::ID temp(0);
    char open, comma, close;
    is >> open >> x.position >> comma >> temp >> close;
    if (!is) return is;

    if (open != '(' ||
        comma != ',' ||
        close != ')')
        throw runtime_error("[operator>>(DNABlock)] Invalid format.");
    x.id = temp;
    return is;
}


bool operator<(const DNABlock& a, const DNABlock& b)
{
    return a.position < b.position ||
           a.position == b.position && a.id < b.id;
}


bool operator>(const DNABlock& a, const DNABlock& b)
{
    return b<a;
}


bool operator==(const DNABlock& a, const DNABlock& b)
{
    return a.position == b.position && a.id == b.id;
}


bool operator!=(const DNABlock& a, const DNABlock& b)
{
    return !(a==b);
}


//
// Chromosome::ID
//


Chromosome::ID::ID(unsigned int encoded)
{
    population = (encoded & 0xf0000000) >> 28;
    individual = (encoded & 0x0fffffc0) >> 6;
    pair = (encoded & 0x3e) >> 1;
    which = (encoded & 1);
}


Chromosome::ID::operator unsigned int() const
{
    const unsigned int maxPopulation = 1<<4;
    const unsigned int maxIndividual = 1<<22;
    const unsigned int maxPair = 1<<5;
    const unsigned int maxWhich = 1<<1;
    
    if (population >= maxPopulation) throw runtime_error("[Chromosome::ID] maxPopulation exceeded.");
    if (individual >= maxIndividual) throw runtime_error("[Chromosome::ID] maxIndividual exceeded.");
    if (pair >= maxPair) throw runtime_error("[Chromosome::ID] maxPair exceeded.");
    if (which >= maxWhich) throw runtime_error("[Chromosome::ID] maxWhich exceeded.");

    return (population << 28) |
           (individual << 6) |
           (pair << 1) | 
           (which);
}


ostream& operator<<(ostream& os, const Chromosome::ID& x)
{
    os << "<" << x.population << "," << x.individual << "," << x.pair << "," << x.which << ">";
    return os;
}


istream& operator>>(istream& is, Chromosome::ID& x)
{
    char open, comma1, comma2, comma3, close;
    is >> open >> x.population >> comma1 >> x.individual >> comma2 >> x.pair >> comma3 >> x.which >> close;
    if (!is) return is;

    if (open != '<' ||
        comma1 != ',' ||
        comma2 != ',' ||
        comma3 != ',' ||
        close != '>')
        throw runtime_error("[operator>>(Chromosome::ID)] Invalid format.");
    
    return is;
}


//
// Chromosome
//


Chromosome::Chromosome(unsigned int id)
{
    blocks_.push_back(DNABlock(0, id));
}


Chromosome::Chromosome(const DNABlocks& blocks)
:   blocks_(blocks)
{}


Chromosome::Chromosome(const Chromosome& x, const Chromosome& y, const vector<unsigned int>& positions)
{
    bool copy_from_x = true; // false == copy from y
    size_t position_previous = 0;

    for (vector<unsigned int>::const_iterator position=positions.begin(); position!=positions.end(); ++position)
    {
        const Chromosome* p = copy_from_x ? &x : &y;
        p->extract_blocks(position_previous, *position, blocks_);

        copy_from_x = !copy_from_x; 
        position_previous = *position;
    }

    const Chromosome* p = copy_from_x ? &x : &y;
    p->extract_blocks(position_previous, numeric_limits<unsigned int>::max(), blocks_);
}


void Chromosome::extract_blocks(unsigned int position_begin, 
                                unsigned int position_end,
                                DNABlocks& result) const
{
    DNABlocks::const_iterator begin = lower_bound(blocks_.begin(), blocks_.end(), DNABlock(position_begin,0));
    DNABlocks::const_iterator end = lower_bound(blocks_.begin(), blocks_.end(), DNABlock(position_end,0));

    if (begin == blocks_.end() || position_begin < begin->position)
    {
        if (begin-1 < blocks_.begin()) throw runtime_error("[Chromosome::extract_blocks()] Blech.");
        result.push_back(DNABlock(position_begin, (begin-1)->id));
    }

    copy(begin, end, back_inserter(result));
}


namespace
{

struct ComparePosition
{
    bool operator()(const DNABlock& a, const DNABlock& b) const {return a.position < b.position;}
};

} // namespace


const DNABlock& Chromosome::find_block(unsigned int position, size_t index_begin)
{
    if (index_begin >= blocks_.size()) throw runtime_error("[Chromosome::find_block()] Bad index_begin.");

    DNABlocks::const_iterator it = upper_bound(blocks_.begin() + index_begin, blocks_.end(), 
                                               DNABlock(position, 0), ComparePosition()); // returns first block with higher position
    return *(it-1);
}


void Chromosome::read(istream& is)
{
    size_t block_count = 0;
    is.read((char*)&block_count, sizeof(size_t));
    if (block_count > 10000) throw runtime_error("[Chromosome::read()] Bad block_count.");
    blocks_.resize(block_count);
    is.read((char*)&blocks_[0], sizeof(DNABlock)*block_count);
}


void Chromosome::write(ostream& os) const
{
    size_t block_count = blocks_.size();
    os.write((const char*)&block_count, sizeof(size_t));
    os.write((const char*)&blocks_[0], sizeof(DNABlock)*block_count);
}


bool operator==(const Chromosome& a, const Chromosome& b)
{
    if (a.blocks().size() != b.blocks().size()) return false;

    for (DNABlocks::const_iterator it=a.blocks().begin(), jt=b.blocks().begin(); it!=a.blocks().end(); ++it,++jt)
        if (*it != *jt) return false;

    return true;
}


bool operator!=(const Chromosome& a, const Chromosome& b)
{
    return !(a==b);
}


ostream& operator<<(ostream& os, const Chromosome& x)
{
    os << "{ ";
    copy(x.blocks().begin(), x.blocks().end(), ostream_iterator<DNABlock>(os, " "));
    os << "}";
    return os;
}


istream& operator>>(istream& is, Chromosome& x)
{
    string buffer;
    getline(is, buffer, '}');
    if (!is) return is;

    istringstream iss(buffer);
    char open;
    iss >> open;
    if (open != '{')
        throw runtime_error("[operator>>(Chromosome)] Invalid format.");

    DNABlocks temp;
    copy(istream_iterator<DNABlock>(iss), istream_iterator<DNABlock>(), back_inserter(temp));
    x = Chromosome(temp);

    return is;
}


