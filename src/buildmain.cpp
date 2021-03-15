#include <iostream>
#include "BOSS/BOSS.hh"
#include "util/cxxopts.hpp"
#include "util/input_reading.hh"
#include "BufBOSS.hh"
#include "Light_merging.hh"
#include "build.hh"
#include "stdafx.h"
#include <iostream>
#include "kmc_file.h"
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

class Kmer_stream_from_files{

    unordered_set<Kmer> kmers;
    unordered_set<Kmer>::iterator it;
    LL k;

public:

    void add_kmers(string& read){
        if((LL)read.size() < k) return;

        Kmer kmer;
        for(LL i = 0; i < k; i++) kmer.appendright(read[i]);
        kmers.insert(kmer);
        for(LL i = k; i < (LL)read.size(); i++){
            kmer.dropleft();
            kmer.appendright(read[i]);
            kmers.insert(kmer);
        }
    }

    Kmer_stream_from_files(const vector<string>& filenames, LL k, bool rc) : k(k){
        
        write_log("Hashing " + to_string(k) + "-mers");
        for(string filename : filenames){
            write_log("Processing file " + filename);
            Sequence_Reader sr(filename, FASTA_MODE);
            while(!sr.done()){
                string read = sr.get_next_query_stream().get_all();
                clean_up_seq(read);
                add_kmers(read);
                if(rc){
                    reverse_complement(read);
                    add_kmers(read);
                }
            }
        }
        write_log(to_string(kmers.size()) + " unique " + to_string(k) + "-mers hashed");
        write_log(to_string(chars_replaced) + " total non-AGCT characters replaced with random nucleotides");
        it = kmers.begin();
    }

    bool done(){
        return it == kmers.end();
    }

    Kmer next(){
        Kmer kmer = *it;
        it++;
        return kmer;
    }
};

class Kmer_stream_from_KMC_DB{

private:

CKMCFile kmer_database;
CKmerAPI kmer_object;

uint32 _kmer_length;
uint32 _mode;
uint32 _counter_size;
uint32 _lut_prefix_length;
uint32 _signature_len;
uint32 _min_count;
uint64 _max_count;
uint64 _total_kmers;

std::string str;
std::string str_revcomp;
bool revcomp_next = false;

public:

    Kmer_stream_from_KMC_DB(string KMC_db_path){
        if (!kmer_database.OpenForListing(KMC_db_path)){
            write_log("Error opening KMC database " + KMC_db_path);
            exit(1);
        }

		kmer_database.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

        kmer_object = CKmerAPI(_kmer_length);
	
    }

    bool done(){
        return !revcomp_next && kmer_database.Eof();
    }

    Kmer next(){
        if(revcomp_next){
            revcomp_next = false;
            return Kmer(str_revcomp);
        }

        float counter_f;
        uint32 counter_i;
		if(_mode){ //quake compatible mode
			kmer_database.ReadNextKmer(kmer_object, counter_f);
		}
		else {			
			kmer_database.ReadNextKmer(kmer_object, counter_i);
		}

        kmer_object.to_string(str);
        str_revcomp = str;
        reverse_complement(str_revcomp);
        if(str != str_revcomp) revcomp_next = true;

        return Kmer(str);

    }

};

template<typename boss_t>
int templated_main(const cxxopts::ParseResult& cli_params){

    string tempdir = cli_params["tempdir"].as<string>();
    string output_dir = cli_params["out"].as<string>();
    string add_file = cli_params["add"].as<string>();
    string add_filelist = cli_params["add-files"].as<string>();
    string KMC_db_path = cli_params["KMC"].as<string>();
    LL k = cli_params["k"].as<LL>();
    bool revcomps = cli_params["revcomp"].as<bool>();
    bool KMC = KMC_db_path != "";

    if(tempdir == ""){
        write_log("Error: temp directory not specified");
        return 1;
    } else temp_file_manager.set_dir(tempdir);

    if(output_dir == ""){
        write_log("Error: output directory not specified");
        return 1;
    }

    if(!KMC && k == 0){
        write_log("Error: k not specified");
    }
    if(!KMC && k >= 32){
        write_log("Error: k must be at most 31.");
        return 1;
    }
    if(KMC && k != 0){
        write_log("Error: When building from KMC, should not give k, because k is defined in the database file");
        return 1;
    }
    if(KMC && revcomps){
        write_log("Error: When building from KMC, should not give -r or --revcomp, because the database file already defines our k-mers.");
        return 1;
    }

    check_dir_exists(output_dir);
    if(add_file != "") check_readable(add_file);
    if(add_filelist != "") check_readable(add_filelist);

    vector<string> addfiles;
    if(add_file != "") addfiles.push_back(add_file);
    if(add_filelist != "")
        for(string line : readlines(add_filelist))
            addfiles.push_back(line);
    if(KMC && (add_file != "" || add_filelist != "")){
        write_log("Error: must give either KMC database or sequence file(s), but not both");
        return 1;
    }

    boss_t boss;

    if(KMC_db_path != ""){
        // Build from KMC
        write_log("Building from KMC database");
        Kmer_stream_from_KMC_DB stream(KMC_db_path);
        BOSS_builder<boss_t, Kmer_stream_from_KMC_DB> builder;
        boss = builder.build(stream);
    } else{
        // Build from fasta
        write_log(to_string(addfiles.size()) + " input files.");
        if(addfiles.size() == 0){
            write_log("No files given");
            return 1;
        }
        Kmer_stream_from_files stream(addfiles, k+1, revcomps);
        BOSS_builder<boss_t, Kmer_stream_from_files> builder;
        boss = builder.build(stream);
    }


    write_log("Built BOSS (k = " + to_string(boss.get_k()) + ") with " + 
            to_string(boss.number_of_nodes()) + " nodes and " + 
            to_string(boss.number_of_edges()) + " edges");
    LL dummy_count = boss.count_dummies();
    write_log("Number of dummy nodes: " + to_string(dummy_count));
    write_log("Number of k-mer nodes: " + to_string(boss.number_of_nodes() - dummy_count));

    write_log("Writing BufBOSS to " + output_dir);
    BufBOSS<boss_t> bufboss(boss);
    bufboss.save_to_disk(output_dir + "/");
    write_log("Done");

    return 0;
}

int main(int argc, char** argv){
    cxxopts::Options options("build","Builds the BOSS data structure");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("k", "Order of the de Bruijn graph. Node are k-mers, edges (k+1)-mers. If building from KMC, don't give this.", cxxopts::value<LL>()->default_value("0"))
      ("o,out", "Output directory.", cxxopts::value<string>()->default_value(""))
      ("a,add", "Path to a fasta-file. Adds all (k+1)-mers of the fasta-file to the index. If building from KMC, don't give this.", cxxopts::value<string>()->default_value(""))
      ("add-files", "Path to a list of fasta-files, one per line. Adds all (k+1)-mers in all the files to the index. If building from KMC, don't give this.", cxxopts::value<string>()->default_value(""))
      ("d,KMC", "Build from KMC database (path to a KMC database). The KMC database consists of two files: xxx.kmc_pre and xxx.kmc_suf. You should give only the xxx part here. The database should be built from canonical k-mers (the default behaviour of KMC)", cxxopts::value<string>()->default_value(""))
      ("r,revcomp", "Include reverse complemented k-mers. If building from KMC, don't give this.", cxxopts::value<bool>()->default_value("false"))
      ("c,rrr", "Use rrr compression on bit vectors.", cxxopts::value<bool>()->default_value("false"))
      ("h,help", "Print instructions.", cxxopts::value<bool>()->default_value("false"))
      ("t,tempdir", "Directory for temporary working space.", cxxopts::value<string>()->default_value(""))
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
