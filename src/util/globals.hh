#pragma once

#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cmath>
#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include <map>
#include <unordered_map>
#include "TempFileManager.hh"
#include <signal.h>
#include "input_reading.hh"
#include "throwing_streams.hh"
#include <chrono>
#include <iomanip>

using namespace std;
using namespace std::chrono;

typedef int64_t LL;
const char read_separator = '$';
Temp_File_Manager temp_file_manager;

long long cur_time_millis(){
	return (std::chrono::duration_cast< milliseconds >(system_clock::now().time_since_epoch())).count();
}

static long long int program_start_millis = cur_time_millis();

double seconds_since_program_start(){
	return (cur_time_millis() - program_start_millis) / 1000.0;
}


string getTimeString(){
    std::time_t result = std::time(NULL);
    string time = std::asctime(std::localtime(&result));
    return time.substr(0,time.size() - 1); // Trim the trailing newline
}

static bool logging_enabled = true;

void enable_logging(){
    logging_enabled = true;
}

void disable_logging(){
    logging_enabled = false;
}

std::mutex write_log_mutex;

void write_log(string message){
    std::lock_guard<std::mutex> lock(write_log_mutex);
    if(logging_enabled){
        std::streamsize default_precision = std::cout.precision();

        std::cerr << 
        std::setprecision(4) << std::fixed <<
        seconds_since_program_start() <<
        std::setprecision(default_precision) << 
        " " << getTimeString() << " " << message << std::endl;
    }
}

map<string,vector<string> > parse_args(int argc, char** argv){
    // Options are argumenta that start with "--". All non-options
    // that come after an option are parameters for that option
    map<string,vector<string> > M; // Option -> list of parameters
    string current_option = "";
    for(LL i = 1; i < argc; i++){
        string S = argv[i];
        if(S.size() >= 2  && S.substr(0,2) == "--"){
            current_option = S;
            M[current_option].resize(0); // Add empty vector for this option.
        } else{
            if(current_option == ""){
                cerr << "Error parsing command line parameters" << endl;
                exit(1);
            }
            M[current_option].push_back(S);
        }
    }
    return M;
}

string figure_out_file_format(string filename){
    for(LL i = filename.size()-1; i >= 0; i--){
        if(filename[i] == '.'){
            string end = filename.substr(i);
            
            if(end == ".fasta") return "fasta";
            if(end == ".fna") return "fasta";
            if(end == ".ffn") return "fasta";
            if(end == ".faa") return "fasta";
            if(end == ".frn") return "fasta";
            if(end == ".fa") return "fasta";

            if(end == ".fastq") return "fastq";
            if(end == ".fq") return "fastq";

            if(end == ".gz") return "gzip";

            throw(runtime_error("Unknown file format: " + filename));
        }
    }
    throw(runtime_error("Unknown file format: " + filename));
    return "unknown";
}

char fix_char(char c){
    char c_new = toupper(c);
    if(c_new != 'A' && c_new != 'C' && c_new != 'G' && c_new != 'T'){
        LL r = rand() % 4;
        if(r == 0) c_new = 'A';
        if(r == 1) c_new = 'C';
        if(r == 2) c_new = 'G';
        if(r == 3) c_new = 'T';
    }
    return c_new;
}

// Returns number of chars replaced
LL fix_alphabet_of_string(string& S){
    LL chars_replaced = 0;
    for(LL i = 0; i < (LL)S.size(); i++){
        char c = S[i];
        char c_new = fix_char(c);
        if(c_new != c){
            S[i] = c_new;
            chars_replaced++;
        }
    }
    return chars_replaced;
}


// Makes a copy of the file and replaces a bad characters. Returns the new filename
// The new file is in fasta format
string fix_alphabet(Sequence_Reader& sr){
    write_log("Making all characters upper case and replacing non-{A,C,G,T} characters with random characeters from {A,C,G,T}");
    //Sequence_Reader fr(fastafile, FASTA_MODE);
    string new_filename = temp_file_manager.get_temp_file_name("seqs-");
    throwing_ofstream out(new_filename);
    LL chars_replaced = 0;
    while(!sr.done()){
        string read = sr.get_next_query_stream().get_all();
        chars_replaced += fix_alphabet_of_string(read);
        out << ">\n" << read << "\n";
    }
    write_log("Replaced " + to_string(chars_replaced) + " characters");
    return new_filename;
}


void sigint_handler(int sig) {
    cerr << "caught signal: " << sig << endl;
    cerr << "Cleaning up temporary files" << endl;
    temp_file_manager.clean_up();
    exit(1);
}

void sigabrt_handler(int sig) {
    cerr << "caught signal: " << sig << endl;
    cerr << "Cleaning up temporary files" << endl;
    temp_file_manager.clean_up();
    cerr << "Aborting" << endl;
    exit(1);
}

auto sigint_register_return_value = signal(SIGINT, sigint_handler); // Set the SIGINT handler
auto sigabrt_register_return_value = signal(SIGABRT, sigabrt_handler); // Set the SIGABRT handler

void set_temp_dir(string temp_dir){
    temp_file_manager.set_dir(temp_dir);
}

string get_temp_file_name(string prefix){
    return temp_file_manager.get_temp_file_name(prefix);
}

vector<string> get_first_and_last_kmers(string fastafile, LL k){
    // todo: this is pretty expensive because this has to read the whole reference data
    Sequence_Reader fr(fastafile, FASTA_MODE);
    vector<string> result;
    while(!fr.done()){
        string ref = fr.get_next_query_stream().get_all();
        if((LL)ref.size() >= k){
            result.push_back(ref.substr(0,k));
            result.push_back(ref.substr(ref.size()-k,k));
        }
    }
    return result;
}



// true if S is colexicographically-smaller than T
bool colex_compare(const string& S, const string& T){
    LL i = 0;
    while(true){
        if(i == (LL)S.size() || i == (LL)T.size()){
            // One of the strings is a suffix of the other. Return the shorter.
            if(S.size() < T.size()) return true;
            else return false;
        }
        if(S[S.size()-1-i] < T[T.size()-1-i]) return true;
        if(S[S.size()-1-i] > T[T.size()-1-i]) return false;
        i++;
    }
}

bool colex_compare_cstrings(const char* x, const char* y){
    LL nx = strlen(x);
    LL ny = strlen(y);
    for(LL i = 0; i < min(nx,ny); i++){
        if(x[nx-1-i] < y[ny-1-i]) return true;
        if(x[nx-1-i] > y[ny-1-i]) return false;
    }

    // All no mismatches -> the shorter string is smaller
    return nx < ny;
};

bool lex_compare(const string& S, const string& T){
    return S < T;
};

bool lex_compare_cstrings(const char* x, const char* y){
    return strcmp(x,y) < 0;
};


template <typename T>
vector<T> parse_tokens(string S){
    vector<T> tokens;
    stringstream ss(S);
    T token;
    while(ss >> token) tokens.push_back(token);
    
    return tokens;
}

// Split by whitespace
vector<string> split(string text){
    std::istringstream iss(text);
    std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                 std::istream_iterator<std::string>());
    return results;
}

// Split by delimiter
vector<string> split(string text, char delimiter){
    assert(text.size() != 0); // If called with empty string we probably have a bug
    vector<LL> I; // Delimiter indices
    I.push_back(-1);
    for(LL i = 0; i < (LL)text.size(); i++){
        if(text[i] == delimiter){
            I.push_back(i);
        }
    }
    I.push_back(text.size());
    vector<string> tokens;
    for(LL i = 0; i < (LL)I.size()-1; i++){
        LL len = I[i+1] - I[i] + 1 - 2;
        tokens.push_back(text.substr(I[i]+1, len));
    }
    
    return tokens;
}

vector<string> split(const char* text, char delimiter){
    return split(string(text), delimiter);
}


// https://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c
void check_dir_exists(string path){
    struct stat info;    
    if( stat( path.c_str(), &info ) != 0 ){
        cerr << "Error: can not access directory " << path << endl;
        exit(1);
    }
    else if( info.st_mode & S_IFDIR ){
        // All good
    }    
    else{
        cerr << "Error: is not a directory: " << path << endl;
        exit(1);
    }
}

void check_readable(string filename){
    throwing_ifstream F(filename); // Throws on failure
}

// Also clears the file
void check_writable(string filename){
    throwing_ofstream F(filename, std::ofstream::out | std::ofstream::app); // Throws on failure
}


template <typename T>
void write_to_file(string path, T& thing){
    throwing_ofstream out(path);
    out << thing << endl;
}


template <typename T>
void read_from_file(string path, T& thing){
    throwing_ifstream input(path);
    input >> thing;
}

vector<string> get_all_lines(string infile){
    vector<string> lines;
    string line;
    throwing_ifstream in(infile);
    while(in.getline(line)){
        lines.push_back(line);
    }
    return lines;
}

vector<char> read_binary_file(string infile){
    throwing_ifstream file(infile, std::ios::binary | std::ios::ate);
    std::streamsize size = file.stream.tellg();
    file.stream.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    if (file.read(buffer.data(), size)){
        return buffer;
    } else{
        cerr << "Error reading file: " << infile << endl;
        assert(false);
    }
    return buffer; // Will never come here but compiler givse a warning if this is not here
}

bool files_are_equal(const std::string& p1, const std::string& p2) {
  //https://stackoverflow.com/questions/6163611/compare-two-files/6163627
    throwing_ifstream f1(p1, std::ifstream::binary|std::ifstream::ate);
    throwing_ifstream f2(p2, std::ifstream::binary|std::ifstream::ate);

    if (f1.stream.tellg() != f2.stream.tellg()) {
      return false; //size mismatch
    }

    //seek back to beginning and use std::equal to compare contents
    f1.stream.seekg(0, std::ifstream::beg);
    f2.stream.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.stream.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.stream.rdbuf()));
}

void check_true(bool condition, string error_message){
    if(!condition){
        throw std::runtime_error(error_message);
    }
}

class Progress_printer{

    public:

    LL n_jobs;
    LL processed;
    LL total_prints;
    LL next_print;

    Progress_printer(LL n_jobs, LL total_prints) : n_jobs(n_jobs), processed(0), total_prints(total_prints), next_print(0) {}

    void job_done(){
        if(next_print == processed){
            LL progress_percent = round(100 * ((double)processed / n_jobs));
            write_log("Progress: " + to_string(progress_percent) + "%");
            next_print += n_jobs / total_prints;
        }
        processed++;
    }

};


set<char> get_alphabet(string S){
    set<char> ans;
    for(char c: S) ans.insert(c);
    return ans;
}

set<string> get_all_distinct_kmers(string S, LL k){
    set<string> kmers;
    for(LL i = 0; i < (LL)(S.size())-k+1; i++){
        kmers.insert(S.substr(i,k));
    }
    return kmers;
}

set<string> get_all_distinct_cyclic_kmers(string S, LL k){
    set<string> kmers;
    for(LL i = 0; i < (LL)S.size(); i++){
        string kmer;
        for(LL j = 0; j < k; j++){
            kmer += S[(i+j) % S.size()];
        }
        kmers.insert(kmer);
    }
    return kmers;
}


set<string> get_all_distinct_cyclic_kmers(vector<string>& A, LL k){
    string concat;
    for(string read : A){
        concat += read_separator + read;
    }
    concat += '\x01'; // bibwt end sentinel

    return get_all_distinct_cyclic_kmers(concat,k);
}

vector<string> get_all_kmers(string& S, LL k){
    vector<string> kmers;
    for(LL i = 0; i < (LL)S.size()-k+1; i++){
        kmers.push_back(S.substr(i,k));
    }
    return kmers;
}

vector<string> all_binary_strings_up_to(int64_t k){ // For testing
    vector<string> ans;
    for(int64_t length = 1; length <= k; length++){
        for(int64_t mask = 0; mask < (1 << length); mask++){
            string s = "";
            for(int64_t i = 0; i < length; i++){
                if(mask & (1 << i)) s += 'A';
                else s += 'C';
            }
            ans.push_back(s);
        }
    }
    return ans;
}

string get_random_dna_string(int64_t length, int64_t alphabet_size){ // For testing
    string s;
    assert(alphabet_size >= 1 && alphabet_size <= 4);
    char alphabet[4] = {'A','T','G','C'};
    for(int64_t i = 0; i < length; i++){
        s.push_back(alphabet[rand() % alphabet_size]);
    }
    return s;
}

string get_random_string(int64_t length, int64_t alphabet_size){ // For testing
    string s;
    for(int64_t i = 0; i < length; i++){
        LL r = rand() % alphabet_size;
        s += 'a' + r;
    }
    return s;
}

vector<string> get_sorted_suffixes(string S){
    vector<string> suffixes;
    for(int64_t i = 0; i < (LL)S.size(); i++){
        suffixes.push_back(S.substr(i));
    }
    sort(suffixes.begin(), suffixes.end());
    return suffixes;
}

template <typename S, typename T>
ostream& operator<<(ostream& os, const unordered_map<S,T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << it->first << ": " << it->second;
    }
    os << "]";
    return os;
}

template <typename S, typename T>
ostream& operator<<(ostream& os, const map<S,T>& v){
    os << "{";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << it->first << ": " << it->second;
    }
    os << "}";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << *it;
    }
    os << "]";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, const set<T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << *it;
    }
    os << "]";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, const multiset<T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << *it;
    }
    os << "]";
    return os;
}

template <typename S, typename T>
ostream& operator<<(ostream& os, const pair<S,T>& x){
    os << "(" << x.first << ", " << x.second << ")";
    return os;
}