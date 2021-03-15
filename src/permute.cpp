#include <iostream>
#include "BOSS/BOSS.hh"
#include "util/cxxopts.hpp"
#include "util/input_reading.hh"
#include "BufBOSS.hh"
#include "Light_merging.hh"
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
    for(char& c : S) c = get_rc(c);
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

    if(index_dir == ""){
        write_log("Error: index directory not given");
        return 1;
    } else {
        write_log("Index directory: " + index_dir);
        check_dir_exists(index_dir);
    }

    write_log("Loading BOSS");
    BufBOSS<boss_t> bufboss;
    bufboss.load_from_disk(index_dir + "/");

    write_log("Benchmarking");
    DynBufferMergeLight merger(bufboss.boss, bufboss.addbuf, bufboss.delbuf, write_log);
    float elapsed = merger.benchmark_one_permutation_iteration();
    write_log("One permutation took " + to_string(elapsed) + " seconds");
    write_log("Done");

    return 0;
    
}

int main(int argc, char** argv){
    cxxopts::Options options("permute", "Benchmarking the node label permutation in BOSS...");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("i,index", "Path to the directory of the index.", cxxopts::value<string>()->default_value(""))
      ("c,rrr", "This option *must* be given if the index was built with rrr compression.", cxxopts::value<bool>()->default_value("false"))
      ("h,help", "Print instructions", cxxopts::value<bool>()->default_value("false"))
    ;

    auto cli_params = options.parse(argc, argv);

    if (cli_params["help"].as<bool>() == true || original_argc == 1){
        cerr << options.help() << endl;
        return 1;
    }

    bool rrr = cli_params["rrr"].as<bool>();
    if(rrr) return templated_main<BOSS<sdsl::rrr_vector<>>>(cli_params);
    return templated_main<BOSS<sdsl::bit_vector>>(cli_params);

}

