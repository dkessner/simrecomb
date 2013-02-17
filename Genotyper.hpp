//
// Genotyper.hpp
//
// Copyright 2013 Darren Kessner
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


#ifndef _GENOTYPER_HPP_
#define _GENOTYPER_HPP_


#include "Population.hpp"
#include "shared_ptr.hpp"
#include <vector>
#include <map>
#include <set>


struct Locus
{
    size_t chromosome_pair_index;
    unsigned int position;

    Locus(size_t _chromosome_pair_index = 0, unsigned int _position = 0) 
    :   chromosome_pair_index(_chromosome_pair_index), position(_position) 
    {}
};


typedef std::set<Locus> Loci;


bool operator<(const Locus& a, const Locus& b);
bool operator==(const Locus& a, const Locus& b);
bool operator!=(const Locus& a, const Locus& b);
std::ostream& operator<<(std::ostream& os, const Locus& locus);


class GenotypeData : public std::vector<char> // {0,1,2} for each individual
{
    public:

    double allele_frequency() const;
};


typedef shared_ptr<GenotypeData> GenotypeDataPtr;
typedef std::map<Locus, GenotypeDataPtr> GenotypeMap;
typedef shared_ptr<GenotypeMap> GenotypeMapPtr;


class SNPIndicator
{
    public:

    // return value in {0, 1}
    virtual unsigned int operator()(unsigned int chromosome_id, const Locus& locus) const = 0;
    virtual ~SNPIndicator() {}
};


typedef shared_ptr<SNPIndicator> SNPIndicatorPtr;


class SNPIndicator_Trivial : public SNPIndicator
{
    public:
    virtual unsigned int operator()(unsigned int chromosome_id, const Locus& locus) const {return 0;}
};


class Genotyper
{
    public:

    // returns genotype value in {0, 1, 2}, for a single organism at a single locus
    unsigned int genotype(const Locus& locus, 
                          const Organism& organism,
                          const SNPIndicator& indicator) const;

    // genotypes population at multiple loci
    GenotypeMapPtr genotype(const Loci& loci, 
                            const Population& population,
                            const SNPIndicator& indicator) const;
};


#endif //  _GENOTYPER_HPP_

