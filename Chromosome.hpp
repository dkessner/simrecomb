//
// Chromosome.hpp
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

#ifndef _CHROMOSOME_HPP_
#define _CHROMOSOME_HPP_


#include <iosfwd>
#include <vector>


struct DNABlock
{
    unsigned int position;       // position in chromosome
    unsigned int id;             // id of original DNA source

    DNABlock(unsigned int _position = 0, unsigned int _id = 0) : position(_position), id(_id) {}
};


typedef std::vector<DNABlock> DNABlocks;


std::ostream& operator<<(std::ostream& os, const DNABlock& x);
std::istream& operator>>(std::istream& is, DNABlock& x);
bool operator<(const DNABlock& a, const DNABlock& b);
bool operator>(const DNABlock& a, const DNABlock& b);
bool operator==(const DNABlock& a, const DNABlock& b);
bool operator!=(const DNABlock& a, const DNABlock& b);


class Chromosome
{
    public:

    // ID structure, encoded as a single int
    struct ID
    {
        unsigned int population; // (4 bits)
        unsigned int individual; // (22 bits)
        unsigned int pair;       // (5 bits) e.g. 0-22 for humans
        unsigned int which;      // (1 bit) 
        
        ID(unsigned int _population,
           unsigned int _individual,
           unsigned int _pair,
           unsigned int _which) 
        :   population(_population),
            individual(_individual),
            pair(_pair),
            which(_which)
        {}

        // decoding from unsigned int
        ID(unsigned int encoded);

        // conversion to unsigned int
        operator unsigned int() const;
    };

    // new chromosome, with single DNABlock
    Chromosome(unsigned int id); 

    // new chromosome with specified DNABlocks, for testing
    Chromosome(const DNABlocks& blocks);

    // new chromosome via recombination
    // note:
    //  - positions must be sorted
    //  - 0 in positions <--> start with y
    Chromosome(const Chromosome& x, const Chromosome& y, const std::vector<unsigned int>& positions); 

    // const access to DNABlocks
    const DNABlocks& blocks() const {return blocks_;}

    // append blocks (from position_begin to position_end) to result
    void extract_blocks(unsigned int position_begin, unsigned int position_end, DNABlocks& result) const;

    private:
    DNABlocks blocks_;
};


typedef std::pair<Chromosome,Chromosome> ChromosomePair;
typedef std::vector<ChromosomePair> ChromosomePairs;


std::ostream& operator<<(std::ostream& os, const Chromosome::ID& x);
std::istream& operator>>(std::istream& is, Chromosome::ID& x);


bool operator==(const Chromosome& a, const Chromosome& b);
bool operator!=(const Chromosome& a, const Chromosome& b);
std::ostream& operator<<(std::ostream& os, const Chromosome& x);
std::istream& operator>>(std::istream& is, Chromosome& x);


#endif //  _CHROMOSOME_HPP_

