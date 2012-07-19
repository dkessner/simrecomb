//
// MSFormat.hpp
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

#ifndef _MSFORMAT_HPP_
#define _MSFORMAT_HPP_


#include <vector>
#include <string>


struct MSFormat
{
    std::vector<double> positions;
    std::vector<std::string> sequences;
    size_t segsites() const {return positions.size();}

    MSFormat(const std::string& filename = "");

    std::string sequence(size_t sequenceIndex, double positionBegin, double positionEnd) const;
};


bool operator==(const MSFormat& a, const MSFormat& b);
bool operator!=(const MSFormat& a, const MSFormat& b);
std::ostream& operator<<(std::ostream& os, const MSFormat& msformat);
std::istream& operator>>(std::istream& is, MSFormat& msformat);


#endif // _MSFORMAT_HPP_

