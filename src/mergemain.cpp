#include <iostream>
#include "BOSS/BOSS.hh"
#include <sys/stat.h>
#include "util/cxxopts.hpp"
#include "util/input_reading.hh"
#include "BufBOSS.hh"
#include "Light_merging.hh"
#include <algorithm>

vector<char> ACGT = {'A','C','G','T'};
LL chars_replaced = 0;
LL min_flush_size = 10000; // Minimum number of k-mers in addition buffer before it is flushed.

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

vector<string> readlines(const string& filename){
    throwing_ifstream in(filename);
    string line;
    vector<string> lines;
    while(in.getline(line)){
        lines.push_back(line);
    }
    return lines;
}

template<typename bufboss_t>
void do_deletions(bufboss_t& bufboss, string fastafile, bool revcomps){
    Sequence_Reader sr(fastafile, FASTA_MODE);
    while(!sr.done()){
        string S = sr.get_next_query_stream().get_all();
        clean_up_seq(S);
        bufboss.delete_edgemers_from_buffer_not_including_dummies(S);
        if(revcomps){
            reverse_complement(S);
            bufboss.delete_edgemers_from_buffer_not_including_dummies(S);
        }
    }
}

template<typename bufboss_t>
void do_additions(bufboss_t& bufboss, string fastafile, bool revcomps, double buffer_fraction){
    Sequence_Reader sr(fastafile, FASTA_MODE);
    while(!sr.done()){
        string S = sr.get_next_query_stream().get_all();
        clean_up_seq(S);
        bufboss.add_edgemers_to_buffer_not_including_dummies(S);
        if(revcomps){
            reverse_complement(S);
            bufboss.add_edgemers_to_buffer_not_including_dummies(S);
        }
        LL n_edges = bufboss.boss.number_of_edges();
        if((LL)bufboss.addbuf.size() > max((LL)(n_edges * buffer_fraction), min_flush_size)){
            write_log("Flushing the buffer");
            bufboss.flush();
        }
    }
}

template<typename boss_t>
int templated_main(const cxxopts::ParseResult& cli_params){

    string input_dir = cli_params["index"].as<string>();
    string output_dir = cli_params["out"].as<string>();
    string add_file = cli_params["add"].as<string>();
    string add_filelist = cli_params["add-files"].as<string>();
    string del_file = cli_params["del"].as<string>();
    string del_filelist = cli_params["del-files"].as<string>();
    LL k = cli_params["k"].as<LL>();
    bool revcomps = cli_params["revcomp"].as<bool>();
    bool end_flush = cli_params["end-flush"].as<bool>();
    double buffer_fraction = cli_params["buffer-fraction"].as<double>();
    bool add_before_del = cli_params["add-before-del"].as<bool>();
    bool count_dummies = cli_params["count-dummies"].as<bool>();

    if(k >= 32){
        write_log("Error: k must be at most 31.");
        return 1;
    }

    if(output_dir == ""){
        write_log("Error: output directory not specified");
        return 1;
    }

    if(input_dir != "") check_dir_exists(input_dir);
    mkdir(output_dir.c_str(),0755);
    check_dir_exists(output_dir);
    if(add_file != "") check_readable(add_file);
    if(add_filelist != "") check_readable(add_filelist);
    if(del_file != "") check_readable(del_file);
    if(del_filelist != "") check_readable(del_filelist);

    vector<string> addfiles;
    if(add_file != "") addfiles.push_back(add_file);
    if(add_filelist != "")
        for(string line : readlines(add_filelist))
            addfiles.push_back(line);    

    vector<string> delfiles;
    if(del_file != "") delfiles.push_back(del_file);
    if(del_filelist != "")
        for(string line : readlines(del_filelist))
            delfiles.push_back(line);    

    write_log(to_string(addfiles.size()) + " input files for addition.");
    write_log(to_string(delfiles.size()) + " input files for deletion.");

    if(addfiles.size() == 0 && delfiles.size() == 0){
        write_log("No additions or deletions given");
        return 1;
    }

    BufBOSS<boss_t> bufboss;
    if(input_dir == ""){
        if(k == 0) {
            write_log("Error: input index not given, so k must be specified for the new index.");
            return 1;
        }
        bufboss = BufBOSS<boss_t>(k);
        write_log("Initialized empty BOSS with k = " + to_string(k));
    } else{
        bufboss = BufBOSS<boss_t>(input_dir + "/");
        k = bufboss.boss.get_k();
        write_log("Loaded BOSS (k = " + to_string(bufboss.boss.get_k()) + ") with " + 
            to_string(bufboss.boss.number_of_nodes()) + " nodes and " +
            to_string(bufboss.boss.number_of_edges()) + " edges");
    }

    if(add_before_del){
        write_log("Loading additions");
        for(string file : addfiles){
            write_log("Adding k-mers of file " + file);
            do_additions(bufboss, file, revcomps, buffer_fraction);
        }
        write_log("Loading additions done");
    }

    write_log("Loading deletions");
    for(string file : delfiles){
        write_log("Deleting k-mers of file " + file);
        do_deletions(bufboss, file, revcomps);
    }
    write_log("Loading deletions done");

    if(!add_before_del){
        write_log("Loading additions");
        for(string file : addfiles){
            write_log("Adding k-mers of file " + file);
            do_additions(bufboss, file, revcomps, buffer_fraction);
        }
        write_log("Loading additions done");
    }

    if(end_flush){
        write_log("Flushing the buffer");
        bufboss.flush();
    }

    write_log(to_string(chars_replaced) + " total non-AGCT characters replaced with random nucleotides");
    
    write_log("Merged BOSS (k = " + to_string(bufboss.boss.get_k()) + ") has " + 
            to_string(bufboss.boss.number_of_nodes()) + " nodes and " + 
            to_string(bufboss.boss.number_of_edges()) + " edges");
    if(count_dummies)
        write_log("Number of dummy nodes: " + to_string(bufboss.boss.count_dummies()));
    write_log("Writing updated bufBOSS to " + output_dir);
    bufboss.save_to_disk(output_dir + "/");
    write_log("Done");

    return 0;
}

int main(int argc, char** argv){
    cxxopts::Options options(argv[0],"Update the bufBOSS data structure");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("i,index", "The directory of the BOSS index. If not given, a new BOSS is built.", cxxopts::value<string>()->default_value(""))
      ("k", "If an input index is not given, a new BOSS is built with this k. Otherwise, this k is ignored.", cxxopts::value<LL>()->default_value("0"))
      ("o,out", "Output directory.", cxxopts::value<string>()->default_value(""))
      ("a,add", "Path to a fasta-file. Adds all (k+1)-mers of the fasta-file to the index.", cxxopts::value<string>()->default_value(""))
      ("add-files", "Path to a list of fasta-files, one per line. Adds all (k+1)-mers in all the files to the index", cxxopts::value<string>()->default_value(""))
      ("add-before-del", "If both additions and deletions are given, the deletions are executed first by default. If you want to execute additions first, give this flag.", cxxopts::value<bool>()->default_value("false"))
      ("d,del", "Path to a fasta-file. Deletes all (k+1)-mers of the fasta-file from the index.", cxxopts::value<string>()->default_value(""))
      ("del-files", "Path to a list of fasta-files, one per line. Deletes all (k+1)-mers in all the files from the index", cxxopts::value<string>()->default_value(""))
      ("r,revcomp", "Include reverse complemented k-mers.", cxxopts::value<bool>()->default_value("false"))
      ("c,rrr", "Use rrr compression on bit vectors.", cxxopts::value<bool>()->default_value("false"))
      ("end-flush", "Flush the buffer at the end before writing to disk.", cxxopts::value<bool>()->default_value("false"))
      ("count-dummies", "Count the number of dummy nodes after the update", cxxopts::value<bool>()->default_value("false"))
      ("b,buffer-fraction", "If this fraction is x and boss has n nodes, then the buffer is flushed when it has max(n*x,"+to_string(min_flush_size)+") k-mers.", cxxopts::value<double>()->default_value("0.01"))
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