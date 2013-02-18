//
// Genotyper.cpp
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


#include "Genotyper.hpp"
#include <iostream>
#include <numeric>


//
// Locus
//


bool operator<(const Locus& a, const Locus& b)
{
    return (a.chromosome_pair_index < b.chromosome_pair_index) ||
           (a.chromosome_pair_index == b.chromosome_pair_index && a.position < b.position);
}


bool operator==(const Locus& a, const Locus& b)
{
    return (a.chromosome_pair_index == b.chromosome_pair_index && a.position == b.position);
}


bool operator!=(const Locus& a, const Locus& b)
{
    return !(a==b);
}


std::ostream& operator<<(std::ostream& os, const Locus& locus)
{
    os << "(" << locus.chromosome_pair_index << "," << locus.position << ")";
    return os;
}


//
// GenotypeData
//


double GenotypeData::allele_frequency() const
{
    double sum = accumulate(begin(), end(), 0.0);
    return sum/size()/2;
}


//
// Genotyper
//


unsigned int Genotyper::genotype(const Locus& locus, 
                                 const Organism& organism,
                                 const SNPIndicator& indicator) const
{
    const ChromosomePair& cp = organism.chromosomePairs()[locus.chromosome_pair_index];
    const DNABlock& block0 = cp.first.find_block(locus.position);
    const DNABlock& block1 = cp.second.find_block(locus.position);
    return indicator(block0.id, locus) + indicator(block1.id, locus);

    // TODO: add support for parallel iteration with a block index hint
}


GenotypeMapPtr Genotyper::genotype(const Loci& loci, 
                                   const Population& population,
                                   const SNPIndicator& indicator) const
{
    GenotypeMapPtr genotype_map(new GenotypeMap);

    for (Loci::const_iterator locus=loci.begin(); locus!=loci.end(); ++locus)
    {
        GenotypeDataPtr genotypes(new GenotypeData);
        genotypes->reserve(population.size());
        const Organisms& organisms = population.organisms();

        for (Organisms::const_iterator organism=organisms.begin(); organism!=organisms.end(); ++organism)
            genotypes->push_back(genotype(*locus, *organism, indicator));

        (*genotype_map)[*locus] = genotypes;
    }

    return genotype_map;
}


