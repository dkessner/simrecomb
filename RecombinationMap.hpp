//
// RecombinationMap.hpp
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

#ifndef _RECOMBINATIONMAP_HPP_
#define _RECOMBINATIONMAP_HPP_


#include <vector>
#include <string>


class RecombinationMap
{
    public:

    // construct with filename "genetic_map_..."
    RecombinationMap(const std::string& filename);

    //
    // HapMap recombination rate 3-column data from files "genetic_map_*":
    //     position COMBINED_rate (cM/Mb) Genetic_Map(cM)
    //
    struct Record
    {
        unsigned int position;
        double combinedRate;
        double geneticMap;

        Record(unsigned int _position = 0,
               double _combinedRate = 0,
               double _geneticMap = 0)
        :   position(_position),
            combinedRate(_combinedRate),
            geneticMap(_geneticMap)
        {}
    };
    
    typedef std::vector<Record> Records;
    Records records() const {return records_;}    

    // return a single random position
    unsigned int random_position();

    // return multiple random positions
    std::vector<unsigned int> random_positions();

    private:
    Records records_;
    std::vector<double> recombinationEventDistribution_;
};


std::istream& operator>>(std::istream& is, RecombinationMap::Record& r);
std::ostream& operator<<(std::ostream& os, const RecombinationMap::Record& r);


#endif // _RECOMBINATIONMAP_HPP_

