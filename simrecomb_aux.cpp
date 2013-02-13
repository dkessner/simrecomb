//
// simrecomb_aux.cpp
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


#include "Population.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>


using namespace std;


int main(int argc, char* argv[])
{
    try
    {
        ostringstream usage;
        usage << "Usage: simrecomb_aux <function> [args]\n";
        usage << endl;
        usage << "Functions:\n";
        usage << "    simrecomb_aux txt2pop filename_in filename_out\n";
        usage << "    simrecomb_aux pop2txt filename_in filename_out\n";
        usage << endl;
        usage << "Darren Kessner\n";
        usage << "John Novembre Lab, UCLA\n";

        string function = argc>1 ? argv[1] : "";

        if (function == "txt2pop")
        {
            if (argc < 4) throw runtime_error(usage.str().c_str());
            string filename_in = argv[2];
            string filename_out = argv[3];

            cout << "reading " << filename_in << endl << flush;
            ifstream is(filename_in.c_str());
            if (!is) throw runtime_error(("[simrecomb_aux] Unable to open file " + filename_in).c_str());
            Population p;
            is >> p;

            cout << "writing " << filename_out << endl << flush;
            ofstream os(filename_out.c_str(), ios::binary);
            p.write(os);
            os.close();
        }
        else if (function == "pop2txt")
        {
            if (argc < 4) throw runtime_error(usage.str().c_str());
            string filename_in = argv[2];
            string filename_out = argv[3];

            cout << "reading " << filename_in << endl << flush;
            ifstream is(filename_in.c_str(), ios::binary);
            if (!is) throw runtime_error(("[simrecomb_aux] Unable to open file " + filename_in).c_str());
            Population p;
            p.read(is);

            cout << "writing " << filename_out << endl << flush;
            ofstream os(filename_out.c_str());
            os << p;
            os.close();
        }
        else
        {
            throw runtime_error(usage.str().c_str());
        }

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


