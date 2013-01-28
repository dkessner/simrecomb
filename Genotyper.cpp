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


unsigned int Genotyper::genotype(const Organism& organism, size_t chromosome_pair_index,
                                 unsigned int position, const SNPIndicator& indicator) const
{
    const ChromosomePair& cp = organism.chromosomePairs()[chromosome_pair_index];
    const DNABlock& block0 = cp.first.find_block(position);
    const DNABlock& block1 = cp.second.find_block(position);
    return indicator(position, block0.id) + indicator(position, block1.id);

    // TODO: add support for parallel iteration with a block index
}

