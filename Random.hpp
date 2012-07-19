//
// Random.hpp
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

#ifndef _RANDOM_HPP_
#define _RANDOM_HPP_


#include "boost/shared_ptr.hpp"


//
// simple random number generator, using Boost.Random;
// method names match Python random module
//
class Random
{
    public:
    
    Random(unsigned int seed = 0);

    // set seed
    void seed(unsigned int value);

    // return random integer N with a <= N <= b
    int randint(int a, int b) const;

    // return random double in [0,1)
    double random() const;

    // return random double in [a,b)
    double uniform(double a, double b) const;

    private:
    class Impl;
    boost::shared_ptr<Impl> impl_;
};


#endif // _RANDOM_HPP_

