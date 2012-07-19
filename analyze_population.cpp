//
// analyze_population.cpp
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
#include <cstring>
#include <fstream>
#include <sstream>
#include <map>
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"


using namespace std;
using boost::shared_ptr;
namespace bfs = boost::filesystem;


struct Config
{
    bfs::path filename;
    bfs::path outputDirectory;
};


void analyzePopulation(const Config& config)
{
    Population::Config dummy;
    Population p(dummy);
    bfs::ifstream is(config.filename);

    cout << "Reading population data.\n";
    is >> p;
    if (p.organisms().empty())
        throw runtime_error("Error reading population data.");

    size_t chromosomePairCount = p.organisms()[0].chromosomePairs().size();

    typedef map<DNABlock, unsigned int> BlockHistogram; // DNABlock -> count
    typedef map<unsigned int, unsigned int> BlockHistogram_2kb; // 2kb position -> count
    typedef map<unsigned int, unsigned int> BlockHistogram_50kb; // 50kb position -> count

    for (size_t i=0; i<chromosomePairCount; ++i)
    {
        BlockHistogram histogram_all;
        BlockHistogram histogram_switch; // recombination between two populations

        BlockHistogram_2kb histogram_all_2kb;
        BlockHistogram_2kb histogram_switch_2kb;
        BlockHistogram_2kb histogram_all_2kb_uniq;
        BlockHistogram_2kb histogram_switch_2kb_uniq;

        BlockHistogram_50kb histogram_all_50kb;
        BlockHistogram_50kb histogram_switch_50kb;
        BlockHistogram_50kb histogram_all_50kb_uniq;
        BlockHistogram_50kb histogram_switch_50kb_uniq;

        for (vector<Organism>::const_iterator it=p.organisms().begin(); it!=p.organisms().end(); ++it)
        {
            vector<const Chromosome*> pairOfChromosomes; // hack to iterate through a pair
            pairOfChromosomes.push_back(&it->chromosomePairs()[i].first);
            pairOfChromosomes.push_back(&it->chromosomePairs()[i].second);
            
            for (vector<const Chromosome*>::const_iterator ch=pairOfChromosomes.begin(); ch!=pairOfChromosomes.end(); ++ch)
            {
                for (DNABlocks::const_iterator block=(*ch)->blocks().begin()+1; block!=(*ch)->blocks().end(); ++block)
                {
                    unsigned int position_2kb = block->position / 2000 * 2000; // round down to nearest 2k
                    unsigned int position_50kb = block->position / 50000 * 50000; // round down to nearest 50k

                    histogram_all[*block]++;
                    histogram_all_2kb[position_2kb]++;
                    histogram_all_50kb[position_50kb]++;
                    if (histogram_all[*block] == 1)
                    {
                        histogram_all_2kb_uniq[position_2kb]++;
                        histogram_all_50kb_uniq[position_50kb]++;
                    }
                    
                    Chromosome::ID current(block->id);
                    Chromosome::ID previous((block-1)->id);
                    if (current.population != previous.population)
                    {
                        histogram_switch[*block]++;
                        histogram_switch_2kb[position_2kb]++;
                        histogram_switch_50kb[position_50kb]++;
                        if (histogram_all[*block] == 1)
                        {
                            histogram_switch_2kb_uniq[position_2kb]++;
                            histogram_switch_50kb_uniq[position_50kb]++;
                        }
                    }
                }
            }
        }

        // compute stats for reporting

        unsigned int totalEvents = 0;
        unsigned int totalSwitches = 0;

        // meta-histograms of events binned by # of occurrences in the population

        map<unsigned int, unsigned int> metahist_all; // # events -> count
        map<unsigned int, unsigned int> metahist_switch; // # switches -> count

        ostringstream filename_block_hist;
        filename_block_hist << "block_hist_" << i << ".txt";
        bfs::ofstream os_block_hist(config.outputDirectory / filename_block_hist.str());

        os_block_hist << "# Chromosome " << i << " histogram\n";
        os_block_hist << "# DNABlock(position,id): events  switches\n";
        for (BlockHistogram::const_iterator it=histogram_all.begin(); it!=histogram_all.end(); ++it)
        {
            unsigned int events = it->second;
            unsigned int switches = histogram_switch[it->first];

            os_block_hist << it->first << ": " << events << " " << switches << endl;

            totalEvents += events;
            totalSwitches += switches;
            
            metahist_all[events]++;
            if (switches) metahist_switch[switches]++;
        }

        unsigned int metahist_all_total = 0;
        unsigned int metahist_switch_total = 0;

        for (map<unsigned int,unsigned int>::const_iterator it=metahist_all.begin(); it!=metahist_all.end(); ++it)
            metahist_all_total += it->second; 

        for (map<unsigned int,unsigned int>::const_iterator it=metahist_switch.begin(); it!=metahist_switch.end(); ++it)
            metahist_switch_total += it->second; 

/*
        // meta-histograms of events binned by # of occurrences per 50kb

        unsigned int totalEvents_50kb = 0;
        unsigned int totalSwitches_50kb = 0;

        map<unsigned int, unsigned int> metahist_all_50kb; // # events -> count
        map<unsigned int, unsigned int> metahist_switch_50kb; // # switches -> count

        ostringstream filename_block_hist_50kb;
        filename_block_hist_50kb << "block_hist_50kb" << i << ".txt";
        bfs::ofstream os_block_hist_50kb(config.outputDirectory / filename_block_hist_50kb.str());

        os_block_hist_50kb << "# Chromosome " << i << " histogram 50kb\n";
        os_block_hist_50kb << "# position  events  switches\n";
        for (BlockHistogram_50kb::const_iterator it=histogram_all_50kb.begin(); it!=histogram_all_50kb.end(); ++it)
        {
            unsigned int events = it->second;
            unsigned int switches = histogram_switch_50kb[it->first];

            os_block_hist_50kb << it->first << " " << events << " " << switches << endl;

            totalEvents_50kb += events;
            totalSwitches_50kb += switches;
            
            metahist_all_50kb[events]++;
            if (switches) metahist_switch_50kb[switches]++;
        }

        unsigned int metahist_all_total_50kb = 0;
        unsigned int metahist_switch_total_50kb = 0;

        for (map<unsigned int,unsigned int>::const_iterator it=metahist_all_50kb.begin(); it!=metahist_all_50kb.end(); ++it)
            metahist_all_total_50kb += it->second; 

        for (map<unsigned int,unsigned int>::const_iterator it=metahist_switch_50kb.begin(); it!=metahist_switch_50kb.end(); ++it)
            metahist_switch_total_50kb += it->second; 
 */        
        // print hist_2kb

        ostringstream filename_hist_2kb;
        filename_hist_2kb << "hist_2kb_" << i << ".txt";
        bfs::ofstream os_hist_2kb(config.outputDirectory / filename_hist_2kb.str());

        os_hist_2kb << "# Chromosome " << i << " histogram 2kb\n";
        os_hist_2kb << "# position  events  switches  events_uniq  switches_uniq\n";
        for (BlockHistogram_2kb::const_iterator it=histogram_all_2kb.begin(); it!=histogram_all_2kb.end(); ++it)
        {
            unsigned int events = it->second;
            unsigned int switches = histogram_switch_2kb[it->first];
            unsigned int events_uniq = histogram_all_2kb_uniq[it->first];
            unsigned int switches_uniq = histogram_switch_2kb_uniq[it->first];
            os_hist_2kb << it->first << " " << events << " " << switches << " "
                        << events_uniq << " " << switches_uniq << endl;
        }

        // print hist_50kb

        ostringstream filename_hist_50kb;
        filename_hist_50kb << "hist_50kb_" << i << ".txt";
        bfs::ofstream os_hist_50kb(config.outputDirectory / filename_hist_50kb.str());

        os_hist_50kb << "# Chromosome " << i << " histogram 50kb\n";
        os_hist_50kb << "# position  events  switches  events_uniq  switches_uniq\n";
        for (BlockHistogram_50kb::const_iterator it=histogram_all_50kb.begin(); it!=histogram_all_50kb.end(); ++it)
        {
            unsigned int events = it->second;
            unsigned int switches = histogram_switch_50kb[it->first];
            unsigned int events_uniq = histogram_all_50kb_uniq[it->first];
            unsigned int switches_uniq = histogram_switch_50kb_uniq[it->first];
            os_hist_50kb << it->first << " " << events << " " << switches << " "
                        << events_uniq << " " << switches_uniq << endl;
        }

        // report stats

        ostringstream filename_events;
        filename_events << "events_" << i << ".txt";
        bfs::ofstream os_events(config.outputDirectory / filename_events.str());

        os_events << "# Chromosome " << i << " recombination events\n";
        os_events << "# total observations: " << totalEvents << endl;
        os_events << "# total unique events: " << metahist_all_total << endl;
        os_events << "# total unique switches: " << metahist_switch_total << endl;
        os_events << "# observationCount uniqueEvents uniqueEventsRelative totalEvents totalEventsRelative uniqueSwitches uniqueSwitchesRelative totalSwitches totalSwitchesRelative\n";

        for (size_t j=1; j<=40; j++)
            os_events << j << " "
                << metahist_all[j] << " " << metahist_all[j]/double(metahist_all_total) << " " 
                << j*metahist_all[j] << " " << j*metahist_all[j]/double(totalEvents) << " " 
                << metahist_switch[j] << " " << metahist_switch[j]/double(metahist_switch_total) << " " 
                << j*metahist_switch[j] << " " << j*metahist_switch[j]/double(totalSwitches) << " " 
                << endl;

/*
        // report stats 50kb

        ostringstream filename_events_50kb;
        filename_events_50kb << "events_50kb_" << i << ".txt";
        bfs::ofstream os_events_50kb(config.outputDirectory / filename_events_50kb.str());

        os_events_50kb << "# Chromosome " << i << " recombination events 50kb\n";
        os_events_50kb << "# total observations: " << totalEvents_50kb << endl;
        os_events_50kb << "# total unique events: " << metahist_all_total_50kb << endl;
        os_events_50kb << "# total unique switches: " << metahist_switch_total_50kb << endl;
        os_events_50kb << "# observationCount uniqueEvents uniqueEventsRelative totalEvents totalEventsRelative uniqueSwitches uniqueSwitchesRelative totalSwitches totalSwitchesRelative\n";

        for (size_t j=1; j<=200; j++)
            os_events_50kb << j << " "
                << metahist_all_50kb[j] << " " << metahist_all_50kb[j]/double(metahist_all_total_50kb) << " " 
                << j*metahist_all_50kb[j] << " " << j*metahist_all_50kb[j]/double(totalEvents_50kb) << " " 
                << metahist_switch_50kb[j] << " " << metahist_switch_50kb[j]/double(metahist_switch_total_50kb) << " " 
                << j*metahist_switch_50kb[j] << " " << j*metahist_switch_50kb[j]/double(totalSwitches_50kb) << " " 
                << endl;
*/
    }
}


Config parseCommandLine(int argc, char* argv[])
{
    if (argc != 3)
    {
        cout << "Usage: analyze_population <filename> <outputdir>\n";
        cout << "\n";
        cout << "Darren Kessner\n";
        cout << "John Novembre Lab, UCLA\n";
        throw runtime_error("");
    }

    Config config;
    config.filename = argv[1];
    config.outputDirectory = argv[2];

    if (!bfs::exists(config.filename))
    {
        ostringstream oss;
        oss << "File not found: " << config.filename;
        throw runtime_error(oss.str());
    }

    if (bfs::exists(config.outputDirectory))
    {
        ostringstream oss;
        oss << "File/directory already exists: " << config.outputDirectory;
        throw runtime_error(oss.str());
    }

    bfs::create_directories(config.outputDirectory);

    return config;
}


int main(int argc, char* argv[])
{
    try
    {
        Config config = parseCommandLine(argc, argv); 
        analyzePopulation(config);
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


