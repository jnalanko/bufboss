#include <iostream>
#include <sys/stat.h>
#include "BOSS/BOSS.hh"
#include "util/cxxopts.hpp"
#include "util/input_reading.hh"
#include "BufBOSS.hh"
#include "Light_merging.hh"
#include <iostream>
#include <unordered_set>
#include <algorithm>

vector<char> ACGT = {'A','C','G','T'};
LL chars_replaced = 0;

char get_rc(char c){
    switch(c){
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: cerr << "Error getting reverse complement from " << c << endl; exit(1);
    }   
}

void reverse_complement(string& S){
    std::reverse(S.begin(), S.end());
    for(char& c : S){
        switch(c){
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            default: break; // Don't change c
        }
    }
}   

bool has_only_ACGT(const string& S){
    for(char c : S){
        if(c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
    }
    return true;
}

void clean_up_seq(string& S){
    for(char& c : S){
        if(c != 'A' && c != 'C' && c != 'G' && c != 'T'){
            c = ACGT[rand()%4];
            chars_replaced++;
        }
    }
}

template<typename boss_t>
int templated_main(const cxxopts::ParseResult& cli_params){

    string index_dir = cli_params["index"].as<string>();
    string outfile = cli_params["out"].as<string>();
    string queryfile = cli_params["query"].as<string>();
    bool revcomps = cli_params["revcomp"].as<bool>();

    if(outfile != ""){
        check_writable(outfile);
    }

    if(index_dir == ""){
        write_log("Error: index directory not given");
        return 1;
    } else {
        write_log("Index directory: " + index_dir);
        check_dir_exists(index_dir);
    }

    if(queryfile == ""){
        write_log("Error: query file note given");
        return 1;
    } else {
        write_log("Query file: " + queryfile);
        check_readable(queryfile);
    }

    throwing_ofstream outfile_stream;

    if(outfile == ""){
        write_log("Outfile not given, outputting to stdout");
    } else {
        write_log("Output file: " + outfile);
        check_writable(outfile);
        outfile_stream.open(outfile);
    }

    BufBOSS<boss_t> bufboss;
    bufboss.load_from_disk(index_dir + "/");
    LL k = bufboss.boss.get_k();

    write_log("Loaded BufBOSS (k = " + to_string(k) + ") with " + 
            to_string(bufboss.boss.number_of_nodes()) + " nodes and " + 
            to_string(bufboss.boss.number_of_edges()) + " edges");
    write_log("Addition buffer has " + to_string(bufboss.addbuf.size()) + " nodemers");

    Sequence_Reader sr(queryfile, FASTA_MODE);
    vector<bool> hits;
    vector<bool> rc_hits;

    LL edgemers_searched = 0;
    LL edgemers_found = 0;

    LL start_time = cur_time_millis();

    while(!sr.done()){
        string read = sr.get_next_query_stream().get_all();

        // Check forward edgemers
        hits = bufboss.search_all_edgemers(read);

        // Check reverse complement edgemers
        if(revcomps){
            reverse_complement(read);
            rc_hits = bufboss.search_all_edgemers(read);
            std::reverse(rc_hits.begin(), rc_hits.end());
            for(LL i = 0; i < (LL)hits.size(); i++)
                hits[i] = hits[i] | rc_hits[i];
        }

        // Record stats
        for(bool b : hits){
            if(b) edgemers_found++;
            edgemers_searched++;
        }

        // Write output
        if(outfile == ""){
            for(bool b : hits) cout << b;
            cout << "\n";
        } else{
            for(bool b : hits) outfile_stream << b;
            outfile_stream << "\n";
        }
    }

    LL elapsed_millis = cur_time_millis() - start_time;
    write_log("Time for all queries: " + to_string((double)elapsed_millis / 1e3) + " seconds");
    write_log(to_string(edgemers_found) + "/" + to_string(edgemers_searched) + " edgemers found (" + 
                        to_string((double)edgemers_found / edgemers_searched * 100) + "%)");
    write_log("Search speed was " + to_string(1e3 * elapsed_millis / edgemers_searched) + " microseconds per edgemer");
    write_log("Done");

    return 0;
}

int main(int argc, char** argv){
    cxxopts::Options options(argv[0],"For every input read R, prints one line L consisting of characters '0' and '1' such L[i] == '1' iff (k+1)-mer R[i..i+k] is found in the index.");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("i,index", "Path to the directory of the index.", cxxopts::value<string>()->default_value(""))
      ("o,out", "Output file. If not given, prints to stdout.", cxxopts::value<string>()->default_value(""))
      ("r,revcomp", "Search reverse-complemented k-mers also.", cxxopts::value<bool>()->default_value("false"))
      ("c,rrr", "This option *must* be given if the index was built with rrr compression.", cxxopts::value<bool>()->default_value("false"))
      ("h,help", "Print instructions.", cxxopts::value<bool>()->default_value("false"))
      ("q,query", "Query FASTA-file", cxxopts::value<string>()->default_value(""))
    ;

    // Todo: save type info in build and merge and use that info to load the right index.

    cxxopts::ParseResult cli_params = options.parse(argc, argv);
    if (cli_params["help"].as<bool>() == true || original_argc == 1){
        cerr << options.help() << endl;
        return 1;
    }

    bool rrr = cli_params["rrr"].as<bool>();
    if(rrr) return templated_main<BOSS<sdsl::rrr_vector<>>>(cli_params);
    return templated_main<BOSS<sdsl::bit_vector>>(cli_params);

}

