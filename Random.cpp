//
// Random.cpp
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
#include "boost/random.hpp"


using namespace std;


struct Random::Impl
{
    boost::mt19937 rng; // generator
    boost::uniform_real<> dist_01; // distribution ~ Uniform(0,1)
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
        random_01; // glues generator to distribution 

    Impl(unsigned int seed)
    :   dist_01(0,1), random_01(rng, dist_01)
    {
        rng.seed(seed);
    }
};


Random::Random(unsigned int seed)
:   impl_(new Impl(seed))
{}


void Random::seed(unsigned int value)
{
    impl_->rng.seed(value);
}


int Random::randint(int a, int b) const
{
/*   
    // not sure if instantiation of these objects takes longer than doing this ourselves
    boost::uniform_int<> dist(a,b); // distribution
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
        vg(impl_->rng, dist); // glues generator to distribution 
    return vg();
*/
    double t = impl_->random_01();
    int result = a + int(t*(b+1-a));
    if (result == b+1) throw runtime_error("[Random::randint()] This isn't happening.");
    return result;
}


double Random::random() const
{
    return impl_->random_01();
}


double Random::uniform(double a, double b) const
{
    return impl_->random_01() * (b-a) + a;
}


