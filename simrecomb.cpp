//
// simrecomb.cpp
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


#include "SimulationController_NeutralAdmixture.hpp"
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
//#include <map>
//#include "boost/filesystem.hpp" // TODO: remove
//#include "boost/filesystem/fstream.hpp" // TODO: remove


using namespace std;
//using boost::shared_ptr;
//namespace bfs = boost::filesystem;


int main(int argc, char* argv[])
{
    try
    {
        ostringstream usage;
        usage << "Usage: simrecomb <subfunction> [args]\n";
        usage << endl;
        usage << "Print example config file to stdout:\n";
        usage << "       simrecomb example\n";
        usage << endl;
        usage << "Run simulation:\n";
        usage << "       simrecomb sim <config_filename> <outputdir>\n";
        usage << endl;
        usage << "Darren Kessner\n";
        usage << "John Novembre Lab, UCLA\n";

        string subfunction = argc>1 ? argv[1] : "";

        if (subfunction == "sim")
        {
            if (argc < 4) throw runtime_error(usage.str().c_str());
                
            string configFilename = argv[2];
            string outputDirectory = argv[3];

            SimulationController::Parameters parameters;
            parameters["popconfig"] = configFilename;
            parameters["outdir"] = outputDirectory;

            SimulationControllerPtr controller(new SimulationController_NeutralAdmixture(parameters));

            controller->initialize();
            controller->run();
        }
        else if (subfunction == "example")
        {
            SimulationControllerPtr controller(new SimulationController_NeutralAdmixture()); // TODO: remove
            controller->example();
        }
        else if (subfunction == "txt2pop")
        {
            if (argc < 4) throw runtime_error(usage.str().c_str());
            string filename_in = argv[2];
            string filename_out = argv[3];

            cout << "reading " << filename_in << endl << flush;
            ifstream is(filename_in.c_str());
            if (!is) throw runtime_error(("[simrecomb] Unable to open file " + filename_in).c_str());
            Population p;
            is >> p;

            cout << "writing " << filename_out << endl << flush;
            ofstream os(filename_out.c_str(), ios::binary);
            p.write(os);
            os.close();
        }
        else if (subfunction == "pop2txt")
        {
            if (argc < 4) throw runtime_error(usage.str().c_str());
            string filename_in = argv[2];
            string filename_out = argv[3];

            cout << "reading " << filename_in << endl << flush;
            ifstream is(filename_in.c_str(), ios::binary);
            if (!is) throw runtime_error(("[simrecomb] Unable to open file " + filename_in).c_str());
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


