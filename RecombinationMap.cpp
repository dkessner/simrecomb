//
// RecombinationMap.cpp
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

#include "RecombinationMap.hpp"
#include "Random.hpp"
#include <iostream>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <cmath>


using namespace std;


namespace {

unsigned int factorial(unsigned int n)
{
    unsigned int result = 1; 
    for (; n>1; n--) result *= n;
    return result;
}

} // namespace


RecombinationMap::RecombinationMap(const string& filename, const Random& random)
:   random_(random)
{
    // read in data file
    ifstream is(filename.c_str());
    string header;
    getline(is, header);
    copy(istream_iterator<RecombinationMap::Record>(is), 
         istream_iterator<RecombinationMap::Record>(), 
         back_inserter(records_));
    if (records_.empty())
        throw runtime_error(("[RecombinationMap] Error reading file " + filename).c_str());

    // calculate Poisson distribution for number of recombination events
    // rate == cumulative geneticMap probability == expected # of events
    double rate = records_.back().geneticMap * .01; // cM * .01 = probability
    double total = 0;
    for (unsigned int i=0; i<10; i++)
    {     
        total += exp(-rate)*pow(rate, double(i))/factorial(i);
        recombinationEventDistribution_.push_back(total);
    }
}


namespace {

struct HasLowerGeneticMap
{
    bool operator()(const RecombinationMap::Record& a, const RecombinationMap::Record& b)
    {
        return a.geneticMap < b.geneticMap;
    }
};

} // namespace


unsigned int RecombinationMap::random_position()
{
    // roll randomly into the distribution, using binary search

    double max = records_.back().geneticMap;
    double roll = random_.uniform(0, max);

    Records::const_iterator it = lower_bound(records_.begin(), records_.end(),
                                             Record(0, 0, roll), HasLowerGeneticMap());

    if (it == records_.begin() || it == records_.end())
        throw runtime_error("[RecombinationMap::random_position()] This isn't happening.");

    // pick a position uniformly between two map positions

    unsigned int range_begin = (it-1)->position;
    unsigned int range_end = it->position - 1;
    unsigned int result = random_.randint(range_begin, range_end);

    return result;
}


vector<unsigned int> RecombinationMap::random_positions()
{
    // random number of events, according to recombinationEventDistribution_

    double roll = random_.random();
    vector<double>::const_iterator it = lower_bound(recombinationEventDistribution_.begin(),
                                                    recombinationEventDistribution_.end(),
                                                    roll);
    size_t count = it - recombinationEventDistribution_.begin();

    // pick random positions

    vector<unsigned int> result;
    for (size_t i=0; i<count; i++)
        result.push_back(random_position());
    return result;
}


istream& operator>>(istream& is, RecombinationMap::Record& r)
{
    is >> r.position >> r.combinedRate >> r.geneticMap;
    return is;
}


ostream& operator<<(ostream& os, const RecombinationMap::Record& r)
{
    os << "(" << r.position << ", " << r.combinedRate << ", " << r.geneticMap << ")";
    return os;
}


