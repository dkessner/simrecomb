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
#include <fstream>
#include <sstream>


using namespace std;


void parse_arg(const string& arg, SimulationController::Parameters& parameters);


void parse_config_file(const string& filename, SimulationController::Parameters& parameters)
{
    cout << "[simrecomb] Parsing config file " << filename << endl; 

    ifstream is(filename.c_str());
    if (!is) throw runtime_error(("[simrecomb] Unable to open file " + filename).c_str());

    while (is)
    {
        string buffer;
        getline(is, buffer);
        if (!buffer.empty() && buffer[0] == '#') continue;
        parse_arg(buffer, parameters); // note: recursion
    }
}


string trim_whitespace(const string& s)
{
    const char* whitespace = " \t\r\n";
    size_t index_begin = s.find_first_not_of(whitespace);
    size_t index_end = s.find_last_not_of(whitespace);
    if (index_begin == string::npos || index_end == string::npos)
        return string();
    return s.substr(index_begin, index_end + 1 - index_begin);
}


void parse_arg(const string& arg, SimulationController::Parameters& parameters)
{
    bool flag = false;

    size_t index_equals = arg.find('=');
    if (index_equals == string::npos)
    {
        index_equals = arg.size();
        flag = true;
    }

    string name = trim_whitespace(arg.substr(0, index_equals));
    if (name.empty()) return;

    size_t index_value = (index_equals < arg.size()) ? index_equals + 1 : index_equals;
    string value = flag ? "1" : trim_whitespace(arg.substr(index_value, arg.size() - index_value));

    //cout << "(name,value): (" << name << "," << value << ")\n";

    if (name == "config")
        parse_config_file(value, parameters);
    else
        parameters[name] = value;
}


void parse_command_line(int argc, char* argv[], string& simname, SimulationController::Parameters& parameters)
{
    ostringstream usage;
    usage << "Usage: simrecomb <simname> [args]\n";
    usage << endl;
    usage << "Available simnames (abbreviation):\n";
    usage << "    neutral_admixture (na)\n";
    usage << endl;
    usage << "Create and run example for <simname>:\n";
    usage << "    simrecomb <simname> example outdir=<dirname>   # creates config files\n";
    usage << "    cd <dirname>\n";
    usage << "    simrecomb <simname> config=config.txt          # runs simulation\n";
    usage << endl;
    usage << "Darren Kessner\n";
    usage << "John Novembre Lab, UCLA\n";

    if (argc < 3)
        throw runtime_error(usage.str().c_str());

    simname = argv[1];
    
    for (int i=2; i<argc; ++i)
        parse_arg(argv[i], parameters);
}


int main(int argc, char* argv[])
{
    try
    {
        string simname;
        SimulationController::Parameters parameters;
        parse_command_line(argc, argv, simname, parameters);

        // instantiate SimulationController

        SimulationControllerPtr controller;
       
        if (simname == "neutral_admixture" || simname == "na")
        {
            controller = SimulationControllerPtr(new SimulationController_NeutralAdmixture(parameters));
        }
        else
        {
            throw runtime_error(("[simrecomb] Unknown simname: "  + simname).c_str());
        }

        if (!controller.get())
            throw runtime_error("[simrecomb] Null SimulationControllerPtr");

        // special handling for "example"

        if (parameters.count("example"))
        {
            if (!parameters.count("outdir"))
                throw runtime_error("[simrecomb] Parameter 'outdir' must be specified for 'example'");

            controller->example(parameters["outdir"]);
            return 0;
        }

        // run the simulation

        controller->initialize();
        controller->run();
        controller->report();

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


