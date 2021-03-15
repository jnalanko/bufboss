#pragma once
#include "BOSS/BOSS.hh"
#include "BOSS/BOSS_construction.hh"
#include <functional>

/**
 * A data record representing the in-edges and out-edges to a k-mer. Used by the merge algorithm.
 */
struct EdgeSet{

    uint16_t data;

    EdgeSet() : data(0) {}

    void set_bit(uint8_t bit, uint8_t value){
        data = (data & ~(((uint16_t)1) << bit)) | (value << bit);
    }

    bool get_bit(uint8_t bit) const{
        return (data >> bit) & ((uint16_t)1);
    }

    bool have_out(char c) const{
        switch(c){
            case 'A': return get_bit(0);
            case 'C': return get_bit(1);
            case 'G': return get_bit(2);
            case 'T': return get_bit(3);
        }
        assert(false); // Invalid character
        return false;
    }

    bool have_in(char c) const{
        switch(c){
            case 'A': return get_bit(4);
            case 'C': return get_bit(5);
            case 'G': return get_bit(6);
            case 'T': return get_bit(7);
        }
        assert(false); // Invalid character
        return false;
    }

    void set_have_out(char c, bool have){
        switch(c){
            case 'A': set_bit(0, have); return;
            case 'C': set_bit(1, have); return;
            case 'G': set_bit(2, have); return;
            case 'T': set_bit(3, have); return;
        }
        assert(false); // Invalid character
    }

    void set_have_in(char c, bool have){
        switch(c){
            case 'A': set_bit(4, have); return;
            case 'C': set_bit(5, have); return;
            case 'G': set_bit(6, have); return;
            case 'T': set_bit(7, have); return;
        }
        assert(false); // Invalid character
    }

    string to_string(){
        stringstream ss;
        ss << "In: {";
        for(char c : string("ACGT")){
            if(have_in(c)) ss << c;
        }
        if(have_indollar()) ss << '$';
        ss << "} Out: {";
        for(char c : string("ACGT")){
            if(have_out(c)) ss << c;
        }
        ss << "}";
        return ss.str();
    }

    static LL size_in_bytes(){
        return sizeof(data);
    }

    void serialize(ostream& out) const{
        out.write((char*)(&data), sizeof(data));
    }

    void load(istream& in){
        in.read((char*)(&data), sizeof(data));
    }

    LL count_out(){
        return (LL)have_out('A') + (LL)have_out('C') + (LL)have_out('G') + (LL)have_out('T');
    }

    LL count_in(){
        return (LL)have_in('A') + (LL)have_in('C') + (LL)have_in('G') + (LL)have_in('T') + (LL)have_indollar();
    }

    bool has_no_edges(){
        return count_in() + count_out() == 0;
    }

    void set_have_indollar(bool have){
        set_bit(8,have);
    }

    bool have_indollar(){
        return get_bit(8);
    }

};

// Wrapper over sdsl::int_vector
// Represents a string of length n with alphabet sigma in space n ceil(log sigma)
class Packed_string{

private:

    sdsl::int_vector<0> v;
    vector<char> alphabet;
    uint8_t char_to_int[256]; // char_to_int[c] = index of character c in the alphabet

public:

    // Helper class that represents a reference to a character in the vector
    struct char_ref{

        const LL idx;
        Packed_string* const S;

        char_ref(LL idx, Packed_string* S) : idx(idx), S(S) {}

        // Cast into a char
        inline operator char() const{
            return S->alphabet[(S->v)[idx]];
        }

        inline char_ref& operator=(char_ref& ref)
        {
            char c = ref;
            this->operator=(c);
            return *this;
        };

        inline char_ref& operator=(char c)
        {
            (S->v)[idx] = S->char_to_int[(uint8_t)c];
            return *this;
        };

        inline bool operator==(const char_ref& other) const{
            return((char)(*this) == (char)other);
        }

        inline bool operator<(const char_ref& other) const{
            return((char)(*this) < (char)other);
        }

        inline bool operator<=(const char_ref& other) const{
            return((char)(*this) <= (char)other);
        }

        inline bool operator>(const char_ref& other) const{
            return((char)(*this) > (char)other);
        }

        inline bool operator>=(const char_ref& other) const{
            return((char)(*this) >= (char)other);
        }

    };

    LL log2_ceiling(LL x){
        
        LL y = 1; // Find the smallest power of 2 that is >= x
        LL log2y = 0;
        while(y < x){
            y *= 2;
            log2y++;
        }
        return log2y;
    }

    Packed_string(const vector<char>& alphabet, LL size, char c) : alphabet(alphabet){
        LL width = max((LL)1,log2_ceiling((LL)alphabet.size()));
        v.width(width);
        v.resize(size);
        for(LL i = 0; i < (LL)alphabet.size(); i++)
            char_to_int[(uint8_t)alphabet[i]] = i;
        for(LL i = 0; i < size; i++){
            v[i] = char_to_int[(uint8_t)c];
        }
    }

    LL size() const {return v.size();}

    // Returns a reference we can assign to
    inline char_ref operator[](LL i){
        char_ref ref(i, this);
        return ref;
    }

    inline char operator[](LL i) const{
        return alphabet[v[i]];
    }

};

// An algorithm that merges a dynamic buffer to a BOSS. To run it, create
// an object of this class and call the merge-function, which will return
// a new BOSS data structure containing both the k-mers in the input BOSS
// and the buffer.
template<typename boss_t>
class DynBufferMergeLight{

template<typename T> friend class DynBufferMergeLightTester;

public:

// Input boss. The result of the algorithm is overwriten to here.
boss_t& boss;

// The addition buffer. Is replaced with an empty unordered_map after the algorithm.
unordered_map<Kmer, EdgeSet>& addbuf_map;

// The delete buffer. Marked in the order of the F-column. Is replaced with an empty vector after the algorithm.
sdsl::bit_vector& delbuf;

private:

    class Null_logger{
        public:
        void operator()(std::string message) {
            (void) message; // Keep compiler happy
            // Do nothing
        }  
    };

    LL n_nodes; // Number of nodes in the input boss
    LL k; // The k of the boss

    // Algorithm internal state
    vector<bool> B1, B2; // Unary coded interval lists. B1 is for boss, B2 is for buf

    // Colex-sorted k-mers and their in/out edge information. This is created from the
    // hash map addition buffer that is given to the constructor.
    vector<pair<Kmer, EdgeSet> > addbuf; 

    // A flag to make sure we don't execute twice on the same object
    bool already_executed = false;

    std::function<void(std::string) > logger = Null_logger();

    // Colex-sort the k-mers in the buffer. WILL CLEAR THE INPUT MAP.
    void get_colex_kmers(){
        for(auto& keyval : addbuf_map) addbuf.push_back(keyval);
        addbuf_map = unordered_map<Kmer, EdgeSet>(); // Free memory
        sort(addbuf.begin(), addbuf.end(), 
            [](const pair<Kmer, EdgeSet>& A, const pair<Kmer, EdgeSet>& B){
                return A.first < B.first;
            }
        );
    }

    // Gets the incoming character of each node
    Packed_string get_first_column(){
        Packed_string col({'\0','A','C','G','T'}, n_nodes, 'A');
        LL n = boss.number_of_nodes() + boss.number_of_edges();
        LL alphabet_size = boss.alphabet_size();
        LL nodes_seen = 0, edges_seen = 0, char_idx = 0;
        for(LL i = 0; i < n; i++){
            if(boss.indegs_at(i) == 1){
                if(i == n-1 || boss.indegs_at(i+1) == 1)
                    col[nodes_seen] = '\0'; // Source node
                else
                    col[nodes_seen] = boss.alphabet_at(char_idx);
                nodes_seen++;
            }
            else{
                edges_seen++;
                if(char_idx < alphabet_size-1 && boss.C_array_at(boss.alphabet_at(char_idx+1)) == edges_seen)
                    char_idx++;
            }
        }
        return col;
    }

    void compute_bitvectors(){

        // Create two full intervals
        B1.resize(boss.number_of_nodes()+1);
        B1[0] = 1;
        B2.resize(addbuf.size()+1);
        B2[0] = 1;

        Packed_string col = get_first_column();
        Packed_string next_col({'\0','A','C','G','T'}, n_nodes, '\0');

        // Take outlabels out the BOSS wavelet tree to speed up access
        Packed_string boss_outlabels({'A','C','G','T'}, boss.number_of_edges(), 'A');// No '\0' needed here
        for(LL i = 0; i < boss.number_of_edges(); i++)
            boss_outlabels[i] = boss.outlabels_at(i);

        for(LL round = 0; round < k; round++){
            logger("Iterating columns. Round " + to_string(round+1) + "/" + to_string(k));
            LL boss_i = 0;
            LL buf_i = 0;
            LL B1_i = 0; LL B2_i = 0;
            vector<bool> B1_new, B2_new;
            while(B1_i < (LL)B1.size() || B2_i < (LL)B2.size()){
                LL len1 = pop_segment(B1, B1_i);
                LL len2 = pop_segment(B2, B2_i);
                LL end1 = boss_i + len1; // One past the end
                LL end2 = buf_i + len2; // One past the end
                if(min(len1,len2) == 0 && max(len1,len2) > 0){
                    push_segment(B1_new, len1); boss_i += len1;
                    push_segment(B2_new, len2); buf_i += len2;
                } else if(min(len1,len2) > 0){
                    // Subdivide: pop runs and pair them up
                    while(boss_i < end1 || buf_i < end2){
                        LL rlen1 = 0; LL rlen2 = 0;
                        if(boss_i == end1){
                            rlen1 = 0;
                            rlen2 = kmer_position_run_length(k-1-round, buf_i, end2);    
                        } else if(buf_i == end2){
                            rlen1 = run_length_from(col, boss_i, end1);
                            rlen2 = 0;
                        } else{
                            char c1 = col[boss_i];
                            char c2 = kmer_col(k-1-round,buf_i);
                            rlen1 = run_length_from(col, boss_i, end1);
                            rlen2 = kmer_position_run_length(k-1-round, buf_i, end2);    
                            if(c1 < c2) rlen2 = 0;
                            if(c1 > c2) rlen1 = 0;
                        }
                        push_segment(B1_new, rlen1); boss_i += rlen1;
                        push_segment(B2_new, rlen2); buf_i += rlen2;
                    }
                }
            }

            permute_boss_labels(col, next_col, boss_outlabels);
            col = next_col;
            B1 = B1_new;
            B2 = B2_new;
        }

    }

    /**
     * Deletions are marked in by the rank of the edge.
     */
    void output_new_boss_with_deletions(){

        const sdsl::bit_vector& deletions = delbuf;

        string new_outlabels;
        vector<bool> new_outdegrees; // Concatenated unary
        vector<bool> new_indegrees; // Concatenated unary

        LL inedge_rank = 0; // Current edge rank in the in-edges of the BOSS.

        auto create_node = [&](){
            new_outdegrees.push_back(1);
            new_indegrees.push_back(1);
        };

        auto output_boss_edges = [&](LL node){
            // This node is in the BOSS but not in the buffer.
            LL l,r;
            std::tie(l,r) = boss.outedge_range(node);
            LL outdegree = r-l+1;
            LL indegree = boss.indegree(node);
            for(LL i = 0; i < outdegree; i++){
                char c = boss.outlabels_at(l+i);
                LL rank = boss.outedge_index_to_wheeler_rank(l+i);
                if(deletions[rank] == 0){ // Not deleted
                    new_outdegrees.push_back(0);
                    new_outlabels.push_back(c);
                }
            }
            for(LL i = 0; i < indegree; i++){
                if(deletions[inedge_rank] == 0){ // Not deleted
                    new_indegrees.push_back(0);
                }
                inedge_rank++;
            }
        };
        
        vector<char> ACGT = {'A','C','G','T'};
        auto output_buffer_edges = [&](LL idx){
            // Add the edges
            for(char c : ACGT){
                if(addbuf[idx].second.have_out(c)){
                    new_outdegrees.push_back(0);
                    new_outlabels.push_back(c);
                }
                if(addbuf[idx].second.have_in(c)){
                    new_indegrees.push_back(0);
                }
            }
            if(addbuf[idx].second.have_indollar()){
                new_indegrees.push_back(0);
            }
        };

        LL B1_i = 0; 
        LL B2_i = 0;

        LL boss_i = 0;
        LL buf_i = 0;
        while(B1_i < (LL)B1.size() || B2_i < (LL)B2.size()){
            LL len1 = pop_segment(B1, B1_i);
            LL len2 = pop_segment(B2, B2_i);
            if(len1 > 0 && len2 > 0){
                // Shared node
                assert(len1 == 1 && len2 == 1);
                create_node();
                output_boss_edges(boss_i++);
                output_buffer_edges(buf_i++);
            } else{
                while(len1--){
                    // Boss only node
                    create_node();
                    output_boss_edges(boss_i++);
                    if(new_indegrees.back() == 1 && new_outdegrees.back() == 1){
                        // No edges were added because they were all marked for deletion
                        // Delete the node. UNLESS it is the source node. We always want to keep that.
                        if(boss_i != 1){
                            new_outdegrees.pop_back();
                            new_indegrees.pop_back();
                        }
                    }
                }
                while(len2--){
                    // Buffer only node
                    create_node();
                    output_buffer_edges(buf_i++);
                }
            }

        }

        delbuf = sdsl::bit_vector(new_outlabels.size(), 0); // Re-initialize
        boss = boss_t(); // Free memory

        // The outlabel sets of a node are now in arbitrary order.
        // Sort them so that we have a canonical representation and we can
        // compare two BOSS-structures bit by bit to see if they are equal
        sort_outlabel_sets(new_outlabels, new_outdegrees);

        // Build the new BOSS support structures
        sdsl::bit_vector sdsl_indegs(new_indegrees.size()), sdsl_outdegs(new_outdegrees.size());
        for(LL i = 0; i < (LL)new_indegrees.size(); i++) sdsl_indegs[i] = new_indegrees[i];
        new_indegrees = vector<bool>(0); // Free memory
        for(LL i = 0; i < (LL)new_outdegrees.size(); i++) sdsl_outdegs[i] = new_outdegrees[i];
        new_outdegrees = vector<bool>(0); // Free memory

        vector<LL> counts(256);
        for(char c : new_outlabels) counts[c]++;
        vector<LL> new_C = char_counts_to_C_array(counts);
        boss = boss_t(new_outlabels, sdsl_indegs, sdsl_outdegs, new_C, k);

    }

    void sort_outlabel_sets(string& outlabels, vector<bool>& outdegrees){
        LL total_num_of_nodes = 0;
        for(LL i = 0; i < (LL)outdegrees.size(); i++) total_num_of_nodes += outdegrees[i];

        LL num_labels_sorted = 0;
        LL outdeg_idx = 0;
        for(LL node = 0; node < total_num_of_nodes; node++){
            LL outdegree = 0;
            outdeg_idx++;
            while(outdeg_idx < (LL)outdegrees.size() && outdegrees[outdeg_idx] == 0){
                outdeg_idx++;
                outdegree++;
            }
            std::sort(outlabels.begin() + num_labels_sorted, outlabels.begin() + num_labels_sorted + outdegree);
            num_labels_sorted += outdegree;
        }
    }

    // Pops a segment like 1 0^s. Returns s. Leaves Bi after the last zero.
    LL pop_segment(const vector<bool>& B, LL& Bi){
        if(Bi == (LL)B.size()) return (LL)0;
        Bi++; // Eat a one-bit
        LL len = 0;
        while(Bi < (LL)B.size() && B[Bi] == 0){
            Bi++; len++;
        }
        return len;
    }

    void push_segment(vector<bool>& B, LL s){
        B.push_back(1);
        while(s--) B.push_back(0);
    }

    // Compute the length of the maximial run of a single element
    // in v starting at position 'start' up to but not including
    // position 'end'
    template<typename vector_t>
    LL run_length_from(const vector_t& v, LL start, LL end){
        if(start >= (LL)v.size()) return (LL)0;
        LL run_end = start+1; // One past the end of the run

        // Repeatedly check if we can add the character at run_end
        while(run_end < min(end, (LL)v.size()) && v[run_end] == v[start]){
            run_end++;
        }
        return run_end-start;
    }

    // Access the character at index pos at buf[idx]
    char kmer_col(LL pos, LL idx){
        // pos = k-1 should give the last character of the k-mer
        LL off = k - addbuf[idx].first.get_k();
        if(pos - off < 0) return '\0';
        return addbuf[idx].first.get(pos - off);
    }

    // Compute the length of the maximial run of a single character
    // at index pos of the k-mers starting from buf[start] up to but
    // not including buf[end].
    LL kmer_position_run_length(LL pos, LL start, LL end){
        if(start >= (LL)addbuf.size()) return (LL)0;
        LL run_end = start+1; // One past the end of the run

        // Repeatedly check if we can add the character at run_end
        while(run_end < min(end, (LL)addbuf.size()) &&
              (kmer_col(pos, run_end) == kmer_col(pos,start))){
                  run_end++;
              }
        return run_end-start;
    }

    // in[i] = inlabel of the p-th predecessor of k-mer i, or null if does not exist
    //         (p = 0) means the inlabel k-mer i
    // result[i] = inlabel of the (p+1)-th predecessor of k-mer i, or null if does not exist
    // outlabels = the outlabels string of the boss. It is given separately to avoid extracting
    // is from the wavelet tree of the boss.
    void permute_boss_labels(Packed_string& in, Packed_string& result, const Packed_string& outlabels){
        for(LL i = 0; i < (LL)result.size(); i++) result[i] = '\0';
        LL idx_outdegrees = 0;
        LL idx_outlabels = 0;
        LL n_nodes = boss.number_of_nodes();
        LL n_edges = boss.number_of_edges();
        map<char, LL> indegs_iterators;
        map<char, LL> indegs_rank1;
        for(LL i = 0; i < boss.alphabet_size(); i++){
            char c = boss.alphabet_at(i);
            indegs_iterators[c] = boss.indegs_select0(boss.C_array_at(c)+1);
            indegs_rank1[c] = boss.indegs_rank1(indegs_iterators[c]);
        }
        
        for(LL node = 0; node < boss.number_of_nodes(); node++){
            
            // Advance in outdegrees bitvector to find out the outdegree
            LL outdegree = 0;
            idx_outdegrees++;
            while(idx_outdegrees < n_nodes + n_edges && boss.outdegs_at(idx_outdegrees) == 0){
                idx_outdegrees++; outdegree++;
            }

            LL outlabels_start = idx_outlabels;
            for(LL i = 0; i < outdegree; i++){
                char c = outlabels[outlabels_start + i];
                // Advance to this edge in indegrees
                while(true){
                    bool bit = boss.indegs_at(indegs_iterators[c]);
                    indegs_iterators[c]++;
                    if(bit == 1) indegs_rank1[c]++;
                    else break;
                }
                LL dest_node = indegs_rank1[c]-1;
                result[dest_node] = in[node];
            }
            idx_outlabels += outdegree;
        }
    }


public:

    typedef std::function<void(std::string) > logger_t;

    /**
     * Sets the references. Running the merge WILL modify them.
     */
    DynBufferMergeLight(
        boss_t& boss, 
        unordered_map<Kmer, EdgeSet>& addbuf, 
        sdsl::bit_vector& delbuf, 
        logger_t logger = Null_logger()) : boss(boss), addbuf_map(addbuf), delbuf(delbuf), n_nodes(boss.number_of_nodes()), k(boss.get_k()) {
            this->logger = logger;
        }

    // The interface to the algorithm. When done, the new boss is in the reference that
    // was given in the constructor. The addition buffer is cleared and the deletion
    // buffer is set to have size equal to number of edgemers in the merged boss. The
    // delete buffer will be filled with all zeroes.
    //
    // Can optinally take a logger callback function. The algorithm will call the logger 
    // with log messages about the progress of the algorithm.
    void merge(){

        assert(!already_executed);
        already_executed = true;

        this->logger("Sorting addition the buffer.");
        get_colex_kmers(); // Clears the addition buffer hash map as a side effect
        logger("Computing the merge plan.");
        compute_bitvectors(); // Compute the merge plan bitvectors B1 and B2
        logger("Executing the merge plan");
        output_new_boss_with_deletions(); // Assigns the new boss to the member variable boss, and re-initializes the delete buffer.
    }

    // Returns number of seconds for running one label permutation
    double benchmark_one_permutation_iteration(){
        logger("Getting first column");
        Packed_string col = get_first_column();

        logger("Initizing buffer for next column");
        Packed_string next_col({'\0','A','C','G','T'}, n_nodes, '\0');

        logger("Getting outlabels out of the rank structure");
        Packed_string boss_outlabels({'A','C','G','T'}, boss.number_of_edges(), 'A'); // No '\0' needed here
        for(LL i = 0; i < boss.number_of_edges(); i++)
            boss_outlabels[i] = boss.outlabels_at(i);

        logger("Starting the timer");
        LL start_time = cur_time_millis();
        permute_boss_labels(col, next_col, boss_outlabels);
        LL end_time = cur_time_millis();
        return (end_time - start_time) / 1000.0;

    }

};