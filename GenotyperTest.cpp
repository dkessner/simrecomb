//
// GenotyperTest.cpp
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
#include "unit.hpp"
#include <iostream>
#include <iterator>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


class SNPIndicator_Test : public SNPIndicator
{
    public:

    virtual unsigned int operator()(unsigned int position, unsigned int id) const
    {
        Chromosome::ID chromosome_id(id);
        return (chromosome_id.population == 0) ? 0 : 1;
    }
};


void test_genotype_easy()
{
    const size_t chromosome_pair_index = 0;
    const unsigned int position = 1000000;
    Locus locus(chromosome_pair_index, position);

    SNPIndicator_Test indicator;
    Chromosome::ID id0(0, 0, 0, 0);
    Chromosome::ID id1(1, 0, 0, 0);
    Genotyper genotyper;

    Organism mom(id0);
    Organism dad(id1);

    unsigned int genotype_mom = genotyper.genotype(locus, mom, indicator);
    unsigned int genotype_dad = genotyper.genotype(locus, dad, indicator);

    if (os_)
    {
        *os_ << "mom:\n" << mom
             << "genotype_mom: " << genotype_mom << endl
             << "dad:\n" << dad
             << "genotype_dad: " << genotype_dad << endl;
    }

    unit_assert(genotype_mom == 0);
    unit_assert(genotype_dad == 2);
}


void test_genotype_harder()
{
    const size_t chromosome_pair_index = 0;
    const unsigned int position = 1000000;
    Locus locus(chromosome_pair_index, position);

    SNPIndicator_Test indicator;
    Chromosome::ID id0(0, 0, 0, 0);
    Chromosome::ID id1(1, 0, 0, 0);
    Genotyper genotyper;

    DNABlocks blocks1;
    blocks1.push_back(DNABlock(0, id0));
    blocks1.push_back(DNABlock(500000, id1));
    blocks1.push_back(DNABlock(1000000, id0)); // SNP 0 in this block
    blocks1.push_back(DNABlock(1500000, id1));
    blocks1.push_back(DNABlock(2000000, id0));

    DNABlocks blocks2;
    blocks2.push_back(DNABlock(0, id1));
    blocks2.push_back(DNABlock(500000, id0));
    blocks2.push_back(DNABlock(900000, id1)); // SNP 1 in this block
    blocks2.push_back(DNABlock(1000001, id0));
    blocks2.push_back(DNABlock(2000000, id1));

    Chromosome chr1(blocks1);
    Chromosome chr2(blocks2);
    
    Organism::Gamete gamete1;
    gamete1.push_back(chr1);

    Organism::Gamete gamete2;
    gamete2.push_back(chr2);

    Organism homo1(gamete1, gamete1);
    Organism homo2(gamete2, gamete2);
    Organism hetero(gamete1, gamete2);

    unsigned int genotype_homo1 = genotyper.genotype(locus, homo1, indicator);
    unsigned int genotype_homo2 = genotyper.genotype(locus, homo2, indicator);
    unsigned int genotype_hetero = genotyper.genotype(locus, hetero, indicator);

    if (os_)
    {
        *os_ << "homo1:\n" << homo1
             << "genotype_homo1: " << genotype_homo1 << endl
             << "homo2:\n" << homo2
             << "genotype_homo2: " << genotype_homo2 << endl
             << "hetero:\n" << hetero
             << "genotype_hetero: " << genotype_hetero << endl;
    }

    unit_assert(genotype_homo1 == 0);
    unit_assert(genotype_homo2 == 2);
    unit_assert(genotype_hetero == 1);

    // test genotype(population)

    Organisms organisms;
    organisms.push_back(homo1);
    organisms.push_back(homo2);
    organisms.push_back(hetero);
    organisms.push_back(homo1);
    organisms.push_back(homo2);
    organisms.push_back(hetero);

    Loci loci;
    loci.push_back(locus);
    
    Population population(organisms);

    GenotypeMapPtr genotype_map = genotyper.genotype(loci, population, indicator);

    GenotypeDataPtr genotypes = genotype_map->at(locus);
    
    if (os_)
    {
        *os_ << "genotypes: " << genotypes->size() << endl;
        copy(genotypes->begin(), genotypes->end(), ostream_iterator<double>(*os_, " "));
        *os_ << endl;
    }

    unit_assert(genotypes->size() == 6);
    unit_assert(genotypes->at(0) == 0 && genotypes->at(3) == 0);
    unit_assert(genotypes->at(1) == 2 && genotypes->at(4) == 2);
    unit_assert(genotypes->at(2) == 1 && genotypes->at(5) == 1);
}


void test()
{
    test_genotype_easy();
    test_genotype_harder();
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


