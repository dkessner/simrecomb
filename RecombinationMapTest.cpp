//
// RecombinationMapTest.cpp
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
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;


void test()
{
    Random random;
    RecombinationMap r("genetic_map_chr21_b36.txt", random);
    unit_assert(r.records().size() == 44250);

    for (int i=0; i<10; i++) 
    {
        vector<unsigned int> positions = r.random_positions();
        if (os_) 
        {
            *os_ << "random_positions(): " << positions.size() << endl;
            copy(positions.begin(), positions.end(), ostream_iterator<unsigned int>(*os_, " "));
            *os_ << endl;
        }
    }
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


