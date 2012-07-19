//
// MSFormat.cpp
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
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <algorithm>


using namespace std;


MSFormat::MSFormat(const string& filename)
{
    if (filename != "")
    {
        ifstream is(filename.c_str());
        if (!is)
            throw runtime_error(("[MSFormat] Unable to open file " + filename).c_str());
        is >> *this;
    }
}


string MSFormat::sequence(size_t sequenceIndex, double positionBegin, double positionEnd) const
{
    if (sequenceIndex >= sequences.size())
        throw runtime_error("[MSFormat::sequence()] Invalid sequence index.");

    vector<double>::const_iterator itBegin = 
        lower_bound(positions.begin(), positions.end(), positionBegin);
    vector<double>::const_iterator itEnd = 
        lower_bound(positions.begin(), positions.end(), positionEnd);

    size_t index = itBegin - positions.begin();
    size_t length = itEnd - itBegin;

    return sequences[sequenceIndex].substr(index, length);
}


bool operator==(const MSFormat& a, const MSFormat& b)
{
    if (a.segsites() != b.segsites()) return false;
    for (vector<double>::const_iterator it=a.positions.begin(), jt=b.positions.begin(); it!=a.positions.end(); ++it, ++jt)
        if (*it != *jt) return false;

    if (a.sequences.size() != b.sequences.size()) return false;
    for (vector<string>::const_iterator it=a.sequences.begin(), jt=b.sequences.begin(); it!=a.sequences.end(); ++it, ++jt)
        if (*it != *jt) return false;

    return true;
}


bool operator!=(const MSFormat& a, const MSFormat& b)
{
    return !(a==b);
}


ostream& operator<<(ostream& os, const MSFormat& msformat)
{
    os << "segsites: " << msformat.segsites() << endl;
    os << "positions: ";
    copy(msformat.positions.begin(), msformat.positions.end(), ostream_iterator<double>(os," "));
    os << endl;
    copy(msformat.sequences.begin(), msformat.sequences.end(), ostream_iterator<string>(os,"\n"));
    return os;
}


istream& operator>>(istream& is, MSFormat& msformat)
{
    msformat.positions.clear();
    msformat.sequences.clear();

    size_t segsites = 0;

    while (is)
    {
        string buffer;
        getline(is, buffer);
        if (!is)
            return is;

        istringstream iss(buffer);
        string first;
        iss >> first;

        if (first == "segsites:")
        {
            iss >> segsites; 
        }
        else if (segsites == 0)
        {
            continue;
        }
        else if (first == "positions:")
        {
            copy(istream_iterator<double>(iss), istream_iterator<double>(), 
                 back_inserter(msformat.positions));
            if (segsites != msformat.segsites())
                throw runtime_error("[operator>>(MSFormat)] Position count doesn't match segsites.");
        }
        else if (!first.empty())
        {
            msformat.sequences.push_back(first);
        }
        else
        {
            break;
        }
    }

    return is;
}


