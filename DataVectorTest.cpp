//
// DataVectorTest.cpp
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


#include "DataVector.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_cdf()
{
    if (os_) *os_ << "test_cdf()\n";

    DataVectorPtr a(new DataVector(5, 1.0));

    if (os_) *os_ << *a << endl;

    DataVectorPtr b = a->cdf();

    if (os_) *os_ << *b << endl;

    for (size_t i=0; i<5; ++i)
        unit_assert(b->at(i) == i+1);
}


void test()
{
    test_cdf();
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


