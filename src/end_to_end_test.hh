#pragma once

#include <string>
#include <unordered_set>
#include "util/input_reading.hh"
#include "BOSS/BOSS.hh"
#include "BufBOSS.hh"
#include "build.hh"
#include "Light_merging.hh"
#include "BOSS/BOSS_tests.hh"
#include <random>
#include <string>
#include <algorithm>

using namespace std;

LL clean_up_seq(string& S){
    vector<char> ACGT = {'A','C','G','T'};
    LL replaced_count = 0;
    for(char& c : S){
        if(c != 'A' && c != 'C' && c != 'G' && c != 'T'){
            c = ACGT[rand()%4];
            replaced_count++;
        }
    }
    return replaced_count;
}

double rand_float(){
    return ((double) rand() / (RAND_MAX));
}

// Picks each string with probability prob
set<string> random_subset(const set<string>& X, double prob){
    set<string> ans;
    for(string S : X) if(rand_float() < prob) ans.insert(S);
    return ans;
}

// Todo: explain what this does
vector<string> segment_string(const string& S, LL edgemer_k){
    vector<bool> bad(S.size());
    LL k = edgemer_k; // Shorthand
    for(LL i = 0; i < (LL)S.size(); i++){
        if(S[i] != 'A' && S[i] != 'C' && S[i] != 'G' && S[i] != 'T'){
            for(LL j = max((LL)0, i-k+1); j <= min((LL)S.size()-1, i+k-1); j++)
                bad[j] = 1;
        }
    }

    vector<string> ans;
    string cur = "";
    for(LL i = 0; i < (LL)S.size(); i++){
        if(!bad[i]) cur += S[i];
        if(bad[i] || i == (LL)S.size() - 1){
            if((LL)cur.size() >= k) ans.push_back(cur);
            cur = "";
        }
    }
    return ans;
}


// Input is the filename of a fasta-file
// Returns a list of sequence such that if a sequence has
// non-ACGT characters, the edgemers containing those are removed
// and the sequence split at the non-ACGT-character.
vector<string> break_up_at_non_ACGT_and_return_all_sequences(string fasta, LL edgemer_k){
    Sequence_Reader sr(fasta, FASTA_MODE);
    vector<string> ans;
    LL n_seqs_unsplitted = 0;
    while(!sr.done()){
        string S = sr.get_next_query_stream().get_all();
        n_seqs_unsplitted++;
        for(string seg : segment_string(S, edgemer_k))
            ans.push_back(seg);
    }
    write_log("Read " + to_string(n_seqs_unsplitted) + " sequences");
    write_log(to_string(ans.size()) + " sequences after splitting at non-ACGT characters");
    return ans;
}

set<string> all_distinct_kmers_of_all_seqs(vector<string>& seqs, LL k){
    set<string> ans;
    for(string S : seqs){
        for(string edgemer : get_all_kmers(S,k)){
            ans.insert(edgemer);
        }
    }
    return ans;
}

void end_to_end_test(string fasta){

    write_log("Running end-to-end test on " + fasta);
    LL k = 30;
    write_log("Edgemers have length " + to_string(k+1));

    write_log("Reading input data");
    vector<string> seqs = break_up_at_non_ACGT_and_return_all_sequences(fasta, k+1);

    LL total_seq_len = 0;
    for(const string& S : seqs) total_seq_len += S.size();

    write_log("Building an initial bufBOSS from the first 50% of total sequence length.");
    vector<Kmer> all_edgemers_packed; // Not all distinct, but *all* edgemers
    for(string S : seqs){
        for(string kmer : get_all_kmers(S,k+1)) all_edgemers_packed.push_back(Kmer(kmer));
    }
    vector<Kmer> half1_edgemers_packed(all_edgemers_packed.begin(), all_edgemers_packed.begin()+all_edgemers_packed.size()/2);
    set<string> true_edgemers; // Distinct edgemers now
    for(Kmer x : half1_edgemers_packed) true_edgemers.insert(x.to_string());

    Kmer_stream_from_iterator_pair<vector<Kmer>::iterator> half1_stream(half1_edgemers_packed.begin(), half1_edgemers_packed.end());
    BOSS_builder<BOSS<sdsl::bit_vector>, Kmer_stream_from_iterator_pair<vector<Kmer>::iterator>> builder;
    BufBOSS bufboss(builder.build(half1_stream));

    write_log("The BOSS of bufboss has " + to_string(bufboss.boss.number_of_nodes()) + " nodes and " + to_string(bufboss.boss.number_of_edges()) + " edges");

    write_log("Getting all our edgemers");
    set<string> our_edgemers = bufboss.boss.get_full_edgemers();

    write_log("Number of true edgemers " + to_string(true_edgemers.size()));
    write_log("Number of our edgemers " + to_string(our_edgemers.size()));
    assert(our_edgemers == true_edgemers);
    write_log("True and our edgemer sets are equal");

    write_log("Adding and deleting all edgemers from all " + to_string(2*k) + "-mers in the whole input data."); 
    write_log("Getting the " + to_string(2*k) + "-mers");

    vector<pair<bool,string> > order; // pairs add/delete, kmer
    for(string S : all_distinct_kmers_of_all_seqs(seqs,2*k)){
        // Let's add it to both BUF+ and BUF- twice so that we get the for sure cases where we add or del the same twice.
        order.push_back({0,S});
        order.push_back({0,S});
        order.push_back({1,S});
        order.push_back({1,S});
    }

    write_log("Shuffling the order of additions and deletions");
    std::shuffle(order.begin(), order.end(), std::default_random_engine(1234));

    write_log("Executing the additions and deletions on the bufboss (no flushing)");
    for(auto op : order){
        if(op.first == 0) bufboss.add_edgemers_to_buffer_not_including_dummies(op.second);
        else bufboss.delete_edgemers_from_buffer_not_including_dummies(op.second);
    }

    write_log("Executing the additions and deletions on the true edgemer set");
    for(auto op : order){
        for(string edgemer : get_all_distinct_kmers(op.second, k+1)){
            if(op.first == 0) true_edgemers.insert(edgemer);
            else true_edgemers.erase(edgemer);
        }
    }
    
    write_log("Checking that every edgemer in the truth is in the bufboss");
    LL n_positive_checks = 0;
    for(string edgemer : true_edgemers){
        assert(bufboss.edgemer_exists(edgemer));
        n_positive_checks++;
    }
    write_log(to_string(n_positive_checks) + " checks passed");

    write_log("Checking that all edgemers of the input data that are not in the truth set anymore are reported not to exist in the bufboss.");
    LL n_negative_checks = 0;
    for(string edgemer : all_distinct_kmers_of_all_seqs(seqs, k+1)){
        if(true_edgemers.count(edgemer) == 0){
            assert(!bufboss.edgemer_exists(edgemer));
            n_negative_checks++;
        }
    }
    write_log(to_string(n_negative_checks) + " checks passed");

    write_log("Check that batch-querying all k-mers of all sequences gives correct answers");

    LL n_batch_checks = 0;
    for(string S : seqs){
        vector<bool> batch_truth;
        for(string edgemer : get_all_kmers(S,k+1)){
            batch_truth.push_back(true_edgemers.count(edgemer));
        }
        vector<bool> our_batch = bufboss.search_all_edgemers(S);
        assert(batch_truth == our_batch);
        n_batch_checks++;
    }
    
    write_log(to_string(n_batch_checks) + " checks passed");

    // Flush and verify
    write_log("Flushing the buffer");
    bufboss.flush();

    write_log("Getting the edgemers of the new BOSS");
    our_edgemers = bufboss.boss.get_full_edgemers();

    write_log("Number of true edgemers " + to_string(true_edgemers.size()));
    write_log("Number of our edgemers " + to_string(our_edgemers.size()));
    assert(our_edgemers == true_edgemers);
    write_log("True and our edgemer sets are equal");

    write_log("Testing counting dummy nodes");
    // Test dummy node counting
    LL dummy_nodes_truth = 0;
    for(LL v = 0; v < (LL)bufboss.boss.number_of_nodes(); v++){
        if((LL)bufboss.boss.get_node_label(v).size() < bufboss.boss.get_k()) dummy_nodes_truth++;
    }
    LL dummy_nodes_our = bufboss.boss.count_dummies();

    write_log("Number of true dummies " + to_string(dummy_nodes_truth));
    write_log("Number of our dummies "  + to_string(dummy_nodes_our));
    assert(dummy_nodes_truth == dummy_nodes_our);
    write_log("Dummy node sets are equal");
    write_log("End-to-end test passes for file " + fasta);

}