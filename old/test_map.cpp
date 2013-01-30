//
// test_map.cpp
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
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <cmath>


using namespace std;


struct GeneticMapRecord
{
    unsigned int position;
    double rate;
    double distance;

    GeneticMapRecord() : position(0), rate(0), distance(0) {}
};


istream& operator>>(istream& is, GeneticMapRecord& r)
{
    is >> r.position >> r.rate >> r.distance;
    return is;
}


ostream& operator<<(ostream& os, const GeneticMapRecord& r)
{
    os << "(" << r.position << ", " << r.rate << ", " << r.distance << ")";
    return os;
}


vector<GeneticMapRecord> read_file(const string& filename)
{
    vector<GeneticMapRecord> result;
    ifstream is(filename.c_str());
    string header;
    getline(is, header);
    copy(istream_iterator<GeneticMapRecord>(is), istream_iterator<GeneticMapRecord>(), back_inserter(result));
    return result;
}


double factorial(double n)
{
    return n==0 ? 1 : n*factorial(n-1);
}


int main()
{
    Random random((int)time(0));
    vector<GeneticMapRecord> v = read_file("genetic_map_chr21_b36.txt");
    //cout << "record count: " << v.size() << endl;
    //copy(v.begin(), v.begin()+10, ostream_iterator<GeneticMapRecord>(cout, "\n"));

    int count = 0;
    double p_total = 0;
    double p0 = 1;
    double f1 = 0, f2 = 0, f3 = 0;

    for (size_t i=1; i<v.size(); ++i)
    {
        double p = .01 * (v[i].distance - v[i-1].distance);

        double roll = random.random();
        bool result = (roll<p);
        //cout << p << " " << roll << " " << boolalpha << result << endl;
        if (result) count++;
        
        p_total += p; 
        p0 *= (1-p);

        double q = p/(1-p);
        f1 += q;
        f2 += q*q;
        f3 += q*q*q;
    }

    cout << "f1: " << f1 << endl;
    cout << "f2: " << f2 << endl;
    cout << "f3: " << f3 << endl;

    cout << "p_total: " << p_total << endl;
    cout << "p(0): " << p0 << endl;
    cout << "p(1): " << p0 * f1 << endl;
    cout << "p(2): " << .5 * p0 * (f1*f1-f2) << endl;
    cout << "p(3): " << p0/6*(f1*f1*f1 - 3*f1*f2 + 2*f3) << endl;

    cout << "Poisson:\n";
    for (size_t i=0; i<7; i++)
        cout << i << " " << exp(-p_total)*pow(p_total, double(i))/factorial(i) << endl;

    cout << "count: " << count << endl;

    return 0;
}



