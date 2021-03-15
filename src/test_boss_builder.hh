#pragma once

#include "util/globals.hh"
#include "BOSS/BOSS_tests.hh"
#include "build.hh"

void test_kmer_stream(vector<string> reads, LL k){
    vector<string> truth;
    for(string S : reads){
        for(LL i = 0; i < (LL)S.size() - k + 1; i++){
            truth.push_back(S.substr(i,k));
        }
    }

    Kmer_stream_in_memory stream(reads, k);
    vector<string> our_edgemers;
    while(!stream.done()){
        our_edgemers.push_back(stream.next().to_string());
    }

    for(LL i = 0; i < (LL)truth.size(); i++){
        assert(our_edgemers[i] == truth[i]);
    }
}

void test_BOSS_builder(){

    cerr << "Testing construction" << endl;

    Boss_Tester<BOSS<sdsl::bit_vector>> tester;

    for(LL k = 1; k <= 12; k++){
        vector<string> reads;
        set<string> true_edgemers;
        for(LL i = 0; i < 100; i++){ // 100*100 = 10k base pairs
            reads.push_back(get_random_dna_string(100,4));
            for(string S : get_all_distinct_kmers(reads.back(), k+1)){
                true_edgemers.insert(S);
            }
        }

        
        //reads = {"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"}; // DEBUG
        // test k-mer stream first
        test_kmer_stream(reads, k+1);

        Kmer_stream_in_memory stream(reads, k+1);
        BOSS_builder<BOSS<sdsl::bit_vector>, Kmer_stream_in_memory> builder;
        BOSS<sdsl::bit_vector> boss = builder.build(stream);

        //cout << tester.BOSS_to_string(boss) << endl;

        set<string> our_edgemers = boss.get_full_edgemers();
        cout << "Our edgemer count " << our_edgemers.size() << endl;
        cout << "True edgemer count " << true_edgemers.size() << endl;
        assert(our_edgemers == true_edgemers);
    }
}