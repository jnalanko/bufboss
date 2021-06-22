#include <iostream>
#include "BOSS/BOSS.hh"
#include "util/cxxopts.hpp"
#include "util/input_reading.hh"
#include "BufBOSS.hh"
#include <iostream>
#include <sys/stat.h>
#include <unordered_set>
#include <algorithm>

template<typename boss_t>
int templated_main(const cxxopts::ParseResult& cli_params){

    string index_dir = cli_params["index"].as<string>();
    string outfile = cli_params["out"].as<string>();

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

    // Flush buffer if there is something to simplify the code that follows
    LL del_count = 0;
    for(LL i = 0; i < (LL)bufboss.delbuf.size(); i++) del_count += bufboss.delbuf[i];
    if(del_count > 0 || bufboss.addbuf.size() > 0)  bufboss.flush();

    write_log("Writing graph");

    // Print graph. Todo: can optimize factor k out of time complexity by detecting and marking all dummy nodes
    for(LL v = 0; v < bufboss.boss.number_of_nodes(); v++){
        if(bufboss.boss.get_node_label(v).size() == bufboss.boss.get_k()){
            for(char c : bufboss.boss.node_outlabels(v)){
                LL u = bufboss.boss.walk(v,c);
                if(outfile == "") cout << v << " " << u << " " << c << "\n";
                else outfile_stream << v << " " << u << " " << c << "\n";
            }
        }
    }



    write_log("Done");

    return 0;
}

int main(int argc, char** argv){
    cxxopts::Options options(argv[0],"Prints the edges of the de Bruijn graph part (no dummy nodes) as integer pairs u,v, one per line");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("i,index", "Path to the directory of the index.", cxxopts::value<string>()->default_value(""))
      ("o,out", "Output file. If not given, prints to stdout.", cxxopts::value<string>()->default_value(""))
      ("c,rrr", "This option *must* be given if the index was built with rrr compression.", cxxopts::value<bool>()->default_value("false"))
      ("h,help", "Print instructions.", cxxopts::value<bool>()->default_value("false"))
    ;

    cxxopts::ParseResult cli_params = options.parse(argc, argv);
    if (cli_params["help"].as<bool>() == true || original_argc == 1){
        cerr << options.help() << endl;
        return 1;
    }

    bool rrr = cli_params["rrr"].as<bool>();
    if(rrr) return templated_main<BOSS<sdsl::rrr_vector<>>>(cli_params);
    return templated_main<BOSS<sdsl::bit_vector>>(cli_params);

}
