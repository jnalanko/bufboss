#pragma once

#include "stdint.h"
#include "Kmer.hh"
#include "Light_merging.hh"
#include <cassert>
#include <bitset>
#include <unordered_map>

using namespace std;

// Dealing with the dummy k-mers is a bit tricky.
// Suppose the edgemers are (k+1)-mers and nodes k-mers.

// When we add an edgemer to the buffer, we add its k-prefix
// and k-suffix, and mark and outedge from the prefix to the suffix.
// 
// If we add an edgemer shorter than k+1, it can be thought as
// a (k+1)-mer that is padded with dollars on the left. Now we
// still add the k-prefix and k-suffix and the edge between them.
// The prefix and suffix may still have dollars. The interpretation
// of a node with dollars is that it is a node with an incoming
// path that terminates at a source node before length (k+1).
//
// The DBG can be seen equivalently as:
//   * A set of nodemers (= k-mers) and a set of edges between them
//   * A set of edgemers (= (k+1)-mers). The nodes are here implicitly as the 
//     prefixes and suffixes of the edgemers. Here is a confusing thing: we
//     consider string that are shorter than (k+1) as edgemers also. It just
//     has an implicit dollar-padding.
//
// Terminology: edgemer is not the same thing as (k+1)-mer. A (k+1)-mer has length
// (k+1). An edgemer can have length \in [1..k+1]. A nodemer can have length [0..k].
// An edgemer is prefixed and suffixed by nodemers. It defines an edge between the
// prefix and the suffix.
//
// Let us fix the order k of the BOSS. Nodemers are up to k and edgemers are up to k+1.
// The BOSS is then fully defined not by the set of nodemers, but the set of 
// edgemers it contains. Not any set of edgemers defines a valid BOSS however.
// It must be the case that: 
// * for every full edgemer S[0..k+1), there is either
//   * A full edgemer T[0..k+1) such that T[1..k+1) = S[0..k), or
//   * A partial edgemer T[0..k) such that T[0..k) = S[0..k)
// * for every partial edgemer S[0..m) with m < k+1, if m != 0, there
//   is a partial edgemer T[0..m-1) = S[0..m-1).
// If these two conditions are satisfied, the BOSS is Wheeler. Proof?
// View the graph in two parts: a tree part and a DBG part. The tree
// part is deterministic, so every node in the tree part has a unique label of length < k,
// so they are sortable and do not interfere with the DBG part. In the DBG part, every node 
// has a unique incoming path suffix of length k, so they are sortable.
//
//

template<typename boss_t = BOSS<sdsl::bit_vector>>
class BufBOSS{

private:
    class Null_logger{
        public:
        void operator()(std::string message) {
            (void) message; // Keep compiler happy
            // Do nothing
        }  
    };

    typedef std::function<void(std::string) > logger_t;
    logger_t logger = Null_logger();

public:

    boss_t boss;
    unordered_map<Kmer, EdgeSet> addbuf;
    sdsl::bit_vector delbuf;

    /**
     *  Constructors an empty bufboss of order 1
     * Can optinally take a logger callback function. The class will call the logger 
     * with log messages about the progress of the algorithm.
     */
    BufBOSS(logger_t logger = Null_logger()) : BufBOSS(1) {
        this->logger = logger;
    }

    /**
     *  Constructors an empty bufboss of order k
     * Can optinally take a logger callback function. The algorithm will call the logger 
     * with log messages about the progress of the algorithm.
     */
    BufBOSS(LL k, logger_t logger = Null_logger()) : boss(k), delbuf(boss.number_of_edges(), 0) {
        this->logger = logger;
    }

    /**
     * Constructs a bufboss with a copy of given boss
     * Can optinally take a logger callback function. The algorithm will call the logger 
     * with log messages about the progress of the algorithm.
     */
    BufBOSS(const boss_t& boss, logger_t logger = Null_logger()) : boss(boss), delbuf(boss.number_of_edges(), 0) {
        this->logger = logger;
    }

    /**
     * Constructs a bufboss by loading a boss from the path prefix
     * Can optinally take a logger callback function. The algorithm will call the logger 
     * with log messages about the progress of the algorithm.
     */
    BufBOSS(string path_prefix, logger_t logger = Null_logger()){
        this->logger = logger;
        load_from_disk(path_prefix);
    }

    void save_to_disk(string path_prefix){
        if(path_prefix.back() != '/') path_prefix += '-';

        // Store boss
        boss.save_to_disk(path_prefix + "boss-");

        // Store delbuf
        sdsl::store_to_file(delbuf, path_prefix + "delbuf");

        // Store addbuf
        throwing_ofstream addbuf_out(path_prefix + "addbuf");
        for(auto keyval : addbuf){
            keyval.first.serialize(addbuf_out.stream);
            keyval.second.serialize(addbuf_out.stream);
        }
    }

    void load_from_disk(string path_prefix){
        if(path_prefix.back() != '/') path_prefix += '-';

        // Load boss
        boss.load_from_disk(path_prefix + "boss-");

        // Load delbuf
        throwing_ifstream delbuf_in(path_prefix + "delbuf");
        delbuf.load(delbuf_in.stream);

        // Load addbuf
        addbuf.clear();
        throwing_ifstream addbuf_in(path_prefix + "addbuf");
        while(true){
            Kmer x; EdgeSet E;
            x.load(addbuf_in.stream);
            E.load(addbuf_in.stream);
            if(addbuf_in.stream.eof()) break;
            addbuf[x] = E;
        }
    }

    /** Add a full (k+1)-mer if it not in the buffer or the boss yet.
     * \return True if addition was succesful, i.e. the edgemer did not already exist.
     */
    bool add_full_edgemer(Kmer edgemer){
        assert(edgemer.get_k() == boss.get_k()+1);
        if(addbuf_has_edgemer(edgemer)) return false; // In addition buffer already
    
        if(boss.edgemer_exists(edgemer.to_string())){ // In boss already
            // Figure out if it is in the delete buffer
            LL rank = boss.edgemer_string_to_edge_wheeler_rank(edgemer.to_string());
            if(delbuf[rank] == 1){
                // Remove from the delete buffer
                delbuf[rank] = 0;
                return true;
            } else return false; // In BOSS but not in delete buffer -> duplicate -> don't add.
        }

        add_full_edgemer_unchecked(edgemer); // Not in boss and not in addition buffer.
        return true;
    }

    /**
     *  Use this only if you know what you are doing
     */
    void add_full_edgemer_unchecked(Kmer edgemer){
        Kmer prefix = edgemer.copy().dropright();
        Kmer suffix = edgemer.copy().dropleft();

        addbuf[prefix].set_have_out(edgemer.get(edgemer.get_k()-1), true);
        addbuf[suffix].set_have_in(edgemer.get(0), true);
    }

    // Use this only if you know what you are doing
    void delete_full_edgemer_from_addition_buffer_unchecked(Kmer edgemer){
        Kmer prefix = edgemer.copy().dropright();
        addbuf[prefix].set_have_out(edgemer.get(edgemer.get_k()-1), false);
        if(addbuf[prefix].has_no_edges())
            addbuf.erase(addbuf.find(prefix));

        Kmer suffix = edgemer.copy().dropleft();
        addbuf[suffix].set_have_in(edgemer.get(0), false);
        if(addbuf[suffix].has_no_edges())
            addbuf.erase(addbuf.find(suffix));
    }

    bool delete_full_edgemer(Kmer edgemer){
        assert(edgemer.get_k() == boss.get_k()+1);
        // Todo optimize: less hash table accesses
        if(addbuf_has_edgemer(edgemer)){
            // In addition buffer
            delete_full_edgemer_from_addition_buffer_unchecked(edgemer);
            return true;
        } else{
            // Not in addition buffer. Find the edge in BOSS
            Kmer prefix = edgemer.copy().dropright();
            LL node = boss.find_kmer(prefix.to_string());
            if(node == -1) return false; // Not in BOSS
            char c = edgemer.get(edgemer.get_k()-1);
            LL nextnode = boss.walk(node, c);
            if(nextnode == -1) return false; // Not in BOSS

            
            // Get the rank of the edge and mark it
            LL edge_idx = boss.outdegs_rank0(boss.outdegs_select1(node+1));
            LL rank = boss.C_array_at(c) + boss.outlabels_rank(edge_idx, c);
            if(delbuf[rank] == 1) return false; // Already in the deletion buffer
            else {
                delbuf[rank] = 1;
                return true;
            }
        }
    }

    /**
     * Returns true if the addition buffer has the edgemer.
     */
    bool addbuf_has_edgemer(Kmer edgemer){
        LL k = edgemer.get_k();
        char c = edgemer.get(k-1);
        edgemer.dropright();
        auto it = addbuf.find(edgemer);
        if(it == addbuf.end()) return false;
        return it->second.have_out(c);
    }

    /**
     * Returns true iff edgemer is in addition buffer or is in BOSS and not in deletion buffer.
     */

    bool edgemer_exists(Kmer edgemer){
        assert(edgemer.get_k() == boss.get_k() + 1);

        if(addbuf_has_edgemer(edgemer)) return true; // In addition buffer

        // Check boss
        LL prev_node = boss.find_kmer(edgemer.to_string().substr(0, boss.get_k()));
        if(prev_node == -1) return false; // Not in BOSS
        LL node = boss.walk(prev_node, edgemer.last());
        if(node == -1) return false; // Not in BOSS

        // Edgemer is in boss. Check deletion buffer.
        LL edge_idx = boss.outdegs_rank0(boss.outdegs_select1(prev_node+1));
        LL rank = boss.C_array_at(edgemer.last()) + boss.outlabels_rank(edge_idx, edgemer.last());
        if(delbuf[rank]) return false; // In BOSS but also in deletion buffer
        return true; // In BOSS and not in deletion buffer
    }

    // Adds all (k+1)-mers of S that are not in the boss to the buffer, not including dummies.
    // k here is the k of the boss.
    void add_edgemers_to_buffer_not_including_dummies(const string& S){
        LL k = boss.get_k();
        // If previous edgemer was found, we don't have to search the next one from scratch.
        // It is enough to try to walk from the previous node
        LL prev_node = -1;
        Kmer edgemer(S.substr(0,k));
        for(LL start = 0; start < (LL)(S.size())-(k+1)+1; start++){
            // Update boss node
            char c = S[start + (k+1) - 1]; // Next character
            LL node = -1;
            if(prev_node == -1) prev_node = boss.find_kmer(S.substr(start,k));
            if(prev_node != -1) node = boss.walk(prev_node, c);

            // Now prev_node is the node of the k-mer starting at S[start] and
            // node is the destination of the (k+1)-mer starting at S[start]. If these nodes
            // don't exist, they are -1.

            if(edgemer.get_k() == k+1) edgemer.dropleft();
            edgemer.appendright(c);
            if(addbuf_has_edgemer(edgemer)){
                // In BUF+ alrady. Nothing to be done
            }
            else if(node != -1){
                // Not in BUF+ but in boss. Remove from BUF- in case it is there.
                LL edge_idx = boss.outdegs_rank0(boss.outdegs_select1(prev_node+1));
                LL rank = boss.C_array_at(c) + boss.outlabels_rank(edge_idx, c);
                delbuf[rank] = 0;
            } else { // Not in BUF+ and not in BOSS -> add to BUF+. Can't be in BUF- because it's not in BOSS.
                add_full_edgemer_unchecked(edgemer);
            } 
            prev_node = node;
        }
    }

    // Deletes all edgemers of S from the buffer.
    void delete_edgemers_from_buffer_not_including_dummies(const string& S){
        LL k = boss.get_k();
        // If previous edgemer was found, we don't have to search the next one from scratch.
        // It is enough to try to walk from the previous node
        LL prev_node = -1;
        Kmer edgemer(S.substr(0,k));
        for(LL start = 0; start < (LL)(S.size())-(k+1)+1; start++){
            // Update boss node
            char c = S[start + (k+1) - 1]; // Next character
            LL node = -1;
            if(prev_node == -1) prev_node = boss.find_kmer(S.substr(start,k));
            if(prev_node != -1) node = boss.walk(prev_node, c);

            // Now prev_node is the node of the k-mer starting at S[start] and
            // node is the destination of the (k+1)-mer starting at S[start]. If these nodes
            // don't exist, they are -1.

            if(edgemer.get_k() == k+1) edgemer.dropleft();
            edgemer.appendright(c);
            if(addbuf_has_edgemer(edgemer)){
                // In BUF+. Remove from BUF+
                delete_full_edgemer_from_addition_buffer_unchecked(edgemer);
            }
            else if(node != -1){
                // Not in BUF+ but in boss
                // Get the rank of the edge and set the delete buffer bit
                LL edge_idx = boss.outdegs_rank0(boss.outdegs_select1(prev_node+1));
                LL rank = boss.C_array_at(c) + boss.outlabels_rank(edge_idx, c);
                delbuf[rank] = 1;
            }
            prev_node = node;
        }
    }

    void flush(){
        LL del_count = 0;
        for(LL i = 0; i < (LL)delbuf.size(); i++) del_count += delbuf[i];
        write_log("Flushing " + to_string(del_count) + " deletions and " + to_string(addbuf.size()) + " additions");
        fill_in_all_dummies();
        DynBufferMergeLight<boss_t> merger(boss, addbuf, delbuf, logger);

        // Run the algorithm. Replaces our boss with a new one. Re-initializes
        // addbuff to empty and delbuf to length equal to the number of edgemers
        // in the new boss.
        merger.merge();
        assert(addbuf.size() == 0);
        assert((LL)delbuf.size() == boss.number_of_edges());
    }

    // Returns a vector v where v[i] = 1 iff edgemer S[i..i+k] contains a character that is not in {A,C,G,T}.
    vector<bool> mark_edgemers_that_have_non_ACGT_chars(const string& S){
        LL k = boss.get_k();
        vector<bool> marks(S.size() - (k+1) + 1);
        LL prev_nonACGT = -1e18;

        for(LL i = 0; i < (LL)S.size(); i++){
            char c = S[i];
            if(c != 'A' && c != 'C' && c != 'G' && c != 'T')
                prev_nonACGT = i;
            if(i >= k){
                if(i - prev_nonACGT < k+1)
                    marks[i-k] = 1;
            }
        }
        return marks;
    }

    vector<bool> search_all_edgemers(const string& S){
        LL k = boss.get_k();
        vector<bool> hits(max((LL)0,(LL)S.size() - (k+1) + 1));
        if(hits.size() == 0) return hits;
        vector<bool> has_non_ACGT = mark_edgemers_that_have_non_ACGT_chars(S);
        // If previous edgemer was found, we don't have to search the next one from scratch.
        // It is enough to try to walk from the previous node
        LL prev_node = -1;
        for(LL start = 0; start < (LL)(S.size())-(k+1)+1; start++){
            if(has_non_ACGT[start]){
                prev_node = -1;
                continue;
            }

            // Update boss node
            char c = S[start + (k+1) - 1]; // Next character
            LL node = -1;
            if(prev_node == -1) prev_node = boss.find_kmer(S.substr(start,k));
            if(prev_node != -1) node = boss.walk(prev_node, c);

            // Now prev_node is the node of the k-mer starting at S[start] and
            // node is the destination of the (k+1)-mer starting at S[start]. If these nodes
            // don't exist, they are -1.

            Kmer edgemer = Kmer(S.substr(start, k+1));
            if(addbuf_has_edgemer(edgemer)){
                hits[start] = 1; // In BUF+.
            }
            else if(node != -1){
                // Not in BUF+ but in boss
                // Get the rank of the edge and check the delete buffer bit
                LL edge_idx = boss.outdegs_rank0(boss.outdegs_select1(prev_node+1));
                LL rank = boss.C_array_at(c) + boss.outlabels_rank(edge_idx, c);
                if(!delbuf[rank]) hits[start] = 1;
            }
            prev_node = node;
        }
        return hits;
    }

private:

    /**
     * Takes in a nodemer (= k-mer).
     * Returns the number of prefixes added 
     */
    LL add_prefixes(const string& S){
        assert((LL)S.size() <= boss.get_k());
        Kmer prefix(0);
        LL node = 0;
        LL addcount = 0;
        for(LL i = 0; i < (LL)S.size(); i++){
            node = boss.walk(node, S[i]);
            prefix.appendright(S[i]);
            if(node == -1){
                if(add_partial_edgemer(prefix)) addcount++;
            }
        }
        return addcount;
    }

    // Filling in dummies. Preconditions:
    // 
    // - Assumes that the buffer has only full nodemers, no dummies. 
    // - Assumes that the edgesets of BOSS, BUF+ and BUF- are all disjoint.
    //
    // If a buffer k-mer does not have an non-deleted incoming edge in BOSS or in BUF+,
    // adds an incoming dollar to the k-mer and adds to the to the buffer 
    // all proper prefixes of the k-mer that are not in BOSS.
    void fill_in_all_dummies(){

        // Make a copy because we don't want to edit the buffer while iterating it.
        // Also make it a vector so it takes less space.
        vector<pair<Kmer, EdgeSet>> old_buf(addbuf.size());
        LL old_buf_idx = 0;
        for(auto& keyval : addbuf) old_buf[old_buf_idx++] = keyval;

        // - For all existing BOSS-nodes that are in BUF+ or at an endpoint of an edge in BUF-: 
        //   find the modified edge set by unioning with BUF+ and subtracting with BUF-. If 
        //   both directions have no edges, then this node will be deleted, so no need for dummies. 
        //   Otherwise, if there are no in-edges, add a dummy chain.
        // - Then, for all *new* nodes in BUF+, if it has no predecessor in BUF+, add a dummy chain.
        //   Note that if a new is new in BUF+, it can not have an in-edge in BOSS, because otherwise
        //   the node would be in the BOSS.
        write_log("Adding dummies");
        write_log("BOSS edgemers: " + to_string(boss.number_of_edges()));
        LL bufplus_edgemers = 0;
        for(auto& keyval : old_buf) bufplus_edgemers += keyval.second.count_out();
        write_log("BUF+ edgemers: " + to_string(bufplus_edgemers));
        LL bufminus_edgemers = 0;
        for(LL i = 0; i < (LL)delbuf.size(); i++) bufminus_edgemers += delbuf[i];
        write_log("BUF- edgemers: " + to_string(bufminus_edgemers));

        LL n_edges_added = 0;
        for(auto& keyval : old_buf){
            Kmer x = keyval.first;
            EdgeSet E = keyval.second;
            assert(x.get_k() == boss.get_k());
            LL node = boss.find_kmer(x.to_string());
            if(node != -1){
                // Old node
                LL out_l, out_r; tie(out_l, out_r) = boss.outedge_range(node);
                LL out_deletes = 0;
                for(LL i = out_l; i <= out_r; i++){
                    LL rank = boss.outedge_index_to_wheeler_rank(i);
                    if(delbuf[rank]) out_deletes++;
                }

                LL in_l, in_r; tie(in_l, in_r) = boss.inedge_range(node);
                LL in_deletes = 0;
                for(LL i = in_l; i <= in_r; i++){
                    if(delbuf[i]) in_deletes++;
                }

                LL new_outdegree = boss.outdegree(node) + E.count_out() - out_deletes;
                LL new_indegree = boss.indegree(node) + E.count_in() - in_deletes;
                if(new_indegree + new_outdegree == 0){
                    // This node will be deleted. No dummies needed
                } else if(new_indegree == 0){
                    n_edges_added += add_prefixes(x.to_string());
                }
            } else{
                // Node in buffer but not BOSS
                // Since it is not in BOSS, there is no incoming edge in BOSS to here.
                if(E.count_in() == 0){
                    // Has no incoming edge in the buffer either.
                    n_edges_added += add_prefixes(x.to_string());
                }
            }
        }

        // Process deletions into nodes that are not in BUF+
        for(LL node = 1; node < boss.number_of_nodes(); node++){
            LL l,r; tie(l,r) = boss.inedge_range(node);
            LL all_incoming_deleted = true;
            for(LL i = l; i <= r; i++) all_incoming_deleted &= delbuf[i];
            if(all_incoming_deleted){
                // Candidate for incoming dummy chain. Check if this node is in BUF+. If yes,
                // it was processed earlier, if no, we need to add a dummy chain, unless the
                // whole node is going to be deleted (i.e. all out-edges are deleted too).
                string label = boss.get_node_label(node);
                if((LL)label.size() == boss.get_k() && addbuf.find(Kmer(label)) == addbuf.end()){
                    LL out_l, out_r; tie(out_l, out_r) = boss.outedge_range(node);
                    LL outdegree = out_r - out_l + 1;
                    LL out_deletes = 0;
                    for(LL i = out_l; i <= out_r; i++){
                        LL rank = boss.outedge_index_to_wheeler_rank(i);
                        if(delbuf[rank]) out_deletes++;
                    }
                    if(outdegree > 0 && out_deletes < outdegree)
                        n_edges_added += add_prefixes(label); // Exists out-edge that is not deleted
                }
            }
        }
        write_log(to_string(n_edges_added) + " dummy edgemers added");
    }    

    /** Add a partial edgemer, i.e. one with length less than k+1,
     * if it not in the buffer or the boss yet.
     * \return True if addition was succesful, i.e. the edgemer did not already exist in the buffer
     */
    bool add_partial_edgemer(Kmer edgemer){
        assert(edgemer.get_k() <= boss.get_k() && edgemer.get_k() >= 1);
        if(addbuf_has_edgemer(edgemer)) return false; // In addition buffer already
        if(boss.edgemer_exists(edgemer.to_string())) return false; // In boss already

        Kmer prefix = edgemer.copy().dropright();
        Kmer suffix = edgemer; // Here we don't drop

        addbuf[prefix].set_have_out(edgemer.get(edgemer.get_k()-1), true);
        addbuf[suffix].set_have_indollar(true);
        return true;
    }

};
