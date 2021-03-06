//
// Organism.hpp
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

#ifndef _ORGANISM_HPP_
#define _ORGANISM_HPP_


#include "Random.hpp"
#include "Chromosome.hpp"
#include "RecombinationMap.hpp"
#include "shared_ptr.hpp"
#include <vector>


//
// RecombinationPositionGenerator interface for generating a list of recombination positions 
//
class RecombinationPositionGenerator
{
    public:
    virtual std::vector<unsigned int> get_positions(size_t index = 0) const = 0;
    virtual ~RecombinationPositionGenerator(){}
};


//
// RecombinationPositionGenerator trivial implementation: returns empty or <0>
// empty == full 1st chromosome, <0> == full 2nd chromosome
//
class RecombinationPositionGenerator_Trivial : public RecombinationPositionGenerator
{
    public:

    RecombinationPositionGenerator_Trivial(const Random& random) : random_(random) {};

    virtual std::vector<unsigned int> get_positions(size_t index) const;

    private:
    const Random& random_;
};


//
// RecombinationPositionGenerator implementation using RecombinationMap
//
class RecombinationPositionGenerator_RecombinationMap : public RecombinationPositionGenerator
{ 
    public:

    RecombinationPositionGenerator_RecombinationMap(const std::vector<std::string>& filenames,
                                                    const Random& random);

    virtual std::vector<unsigned int> get_positions(size_t index) const;

    private:
    const Random& random_;
    std::vector< shared_ptr<RecombinationMap> > recombinationMaps_;
};


class Organism
{
    public:

    static shared_ptr<RecombinationPositionGenerator> recombinationPositionGenerator_;

    typedef std::vector<Chromosome> Gamete;

    Organism(unsigned int id = 0, size_t chromosomeCount = 1);
    Organism(const Gamete& g1, const Gamete& g2);
    Organism(const Organism& mom, const Organism& dad);

    const ChromosomePairs& chromosomePairs() const {return chromosomePairs_;}

    Gamete create_gamete() const;

    // binary read/write
    void read(std::istream& is);
    void write(std::ostream& os) const;

    private:
    ChromosomePairs chromosomePairs_;
};


typedef std::vector<Organism> Organisms;


bool operator==(const Organism& a, const Organism& b);
bool operator!=(const Organism& a, const Organism& b);
std::ostream& operator<<(std::ostream& os, const Organism& o);
std::istream& operator>>(std::istream& is, Organism& o);


#endif // _ORGANISM_HPP_

