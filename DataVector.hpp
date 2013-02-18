//
// DataVector.hpp
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


#ifndef _DATAVECTOR_HPP_
#define _DATAVECTOR_HPP_


#include "shared_ptr.hpp"
#include <iosfwd>
#include <vector>


class DataVector;
typedef shared_ptr<DataVector> DataVectorPtr;
typedef std::vector<DataVectorPtr> DataVectorPtrs;


class DataVector : public std::vector<double>
{
    public:

    DataVector(size_t n = 0, double value = 0)
    :   std::vector<double>(n, value)
    {}

    // convenience functions

    double mean() const;
    DataVectorPtr cdf() const; // note: allocates a new DataVector

    // TODO add functions as needed:  sum, sum_squares, variance
    // inner product -> sum_squares -> variance
};


std::ostream& operator<<(std::ostream& os, const DataVector& v);


#endif //  _DATAVECTOR_HPP_

