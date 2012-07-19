//
// RandomTest.cpp
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

#include "Random.hpp"
#include "unit.hpp"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cstring>
#include <ctime>


using namespace std;


ostream* os_ = 0;


void test()
{
    Random random(static_cast<unsigned int>(std::time(0)));

    if (os_) *os_ << "random():\n";
    for (int i=0; i<10; i++)
        if (os_) *os_ << random.random() << endl;
    if (os_) *os_ << endl;

    if (os_) *os_ << "randint(1,10):\n";
    for (int i=0; i<10; i++)
        if (os_) *os_ << random.randint(1,10) << endl;
    if (os_) *os_ << endl;

    if (os_) *os_ << "uniform(5,7):\n";
    for (int i=0; i<10; i++)
        if (os_) *os_ << random.uniform(5,7) << endl;
    if (os_) *os_ << endl;
}


void test_seed()
{
    if (os_) *os_ << "test_seed()\n";

    Random random(420);

    if (os_) *os_ << "first pass:\n";
    vector<int> v;
    for (int i=0; i<10; i++)
    {
        v.push_back(random.randint(1,10));
        if (os_) *os_ << v.back() << endl;
    }
    if (os_) *os_ << endl;
    
    random.seed(420);

    if (os_) *os_ << "second pass:\n";
    vector<int> w;
    for (int i=0; i<10; i++)
    {
        w.push_back(random.randint(1,10));
        if (os_) *os_ << w.back() << endl;
    }
    if (os_) *os_ << endl;

    for (int i=0; i<10; i++)
        unit_assert(v[i] == w[i]);
}


int main(int argc, char* argv[])
{
    try
    {
        if (argc>1 && !strcmp(argv[1],"-v")) os_ = &cout;
        test();
        test_seed();
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


