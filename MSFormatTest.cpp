//
// MSFormatTest.cpp
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

#include "MSFormat.hpp"
#include "unit.hpp"
#include <iostream>
#include <fstream>
#include <cstring>


using namespace std;


ostream* os_ = 0;


void test0(const string& filename)
{
    ifstream is(filename.c_str());
    unit_assert(is);

    MSFormat ms;
    is >> ms;
    if (os_) *os_ << ms << endl;
    unit_assert(ms.segsites() == 3);
    unit_assert(ms.sequences.size() == 5);
    unit_assert(ms.positions[0] == .3203);
    unit_assert(ms.positions[1] == .6540);
    unit_assert(ms.positions[2] == .8722);

    is >> ms;
    if (os_) *os_ << ms << endl;
    unit_assert(ms.segsites() == 7);
    unit_assert(ms.sequences.size() == 5);
    unit_assert(ms.sequences[4] == "1101100");
}


void test1(const string& filename)
{
    MSFormat ms(filename);
    if (os_) *os_ << ms << endl;
    unit_assert(ms.segsites() == 5359);
    unit_assert(ms.sequences.size() == 100);
}


void testCycle()
{
    if (os_) *os_ << "testCycle()\n";

    MSFormat ms;
    ms.positions.push_back(.25);
    ms.positions.push_back(.5);
    ms.positions.push_back(.75);
    ms.sequences.push_back("000");
    ms.sequences.push_back("001");
    ms.sequences.push_back("010");
    ms.sequences.push_back("011");

    ostringstream oss;
    oss << ms;
    if (os_) *os_ << "ms:\n" << ms << endl;

    MSFormat test;
    unit_assert(test != ms);

    istringstream iss(oss.str());
    iss >> test;
    if (os_) *os_ << "test:\n" << test << endl;
    unit_assert(test == ms);
}


void testSequence()
{
    if (os_) *os_ << "testCycle()\n";

    MSFormat ms;
    ms.positions.push_back(.25);
    ms.positions.push_back(.5);
    ms.positions.push_back(.75);
    ms.sequences.push_back("abc");
    ms.sequences.push_back("def");
    ms.sequences.push_back("ghi");
    ms.sequences.push_back("jkl");

    unit_assert(ms.sequence(0, 0, .3) == "a");
    unit_assert(ms.sequence(0, 0, .5) == "a");
    unit_assert(ms.sequence(0, .5, .75) == "b");
    unit_assert(ms.sequence(0, .25, .75) == "ab");
    unit_assert(ms.sequence(0, .15, .76) == "abc");
    unit_assert(ms.sequence(0, .5, .76) == "bc");
    unit_assert(ms.sequence(1, .15, .85) == "def");
}


int main(int argc, char* argv[])
{
    try
    {
        vector<string> filenames;
        for (int i=1; i<argc; i++)
        {            
            if (!strcmp(argv[i],"-v")) os_ = &cout;
            else filenames.push_back(argv[i]);
        }
        unit_assert(filenames.size() == 2);
    
        if (filenames.size() == 2)
        {
            test0(filenames[0]);
            test1(filenames[1]);
        }

        testCycle();
        testSequence();

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


