#pragma once

#include "Kmer.hh"
#include "util/globals.hh"

void Kmer_basic_tests(){
    string S = get_random_dna_string(32, 4);
    Kmer kmer(32);

    for(LL i = 0; i < 32; i++){
        kmer.set(i, S[i]);
    }

    for(LL i = 0; i < 32; i++){
        assert(kmer.get(i) == S[i]);
    }

    Kmer left = kmer.copy().dropleft();
    assert(left.get_k() == 31);
    for(LL i = 0; i < 31; i++){
        assert(left.get(i) == S[i+1]);
    }

    Kmer right = kmer.copy().dropright();
    assert(right.get_k() == 31);
    for(LL i = 0; i < 31; i++){
        assert(right.get(i) == S[i]);
    }
}

void Kmer_colex_tests(){
    vector<string> strings;
    vector<Kmer> kmers;

    int n = 200;

    // Generate data
    for(int i = 0; i < n; i++){
        strings.push_back(get_random_dna_string(rand() % 4, 4)); // Random length, alphabet size 4
        kmers.push_back(strings.back());
    }

    // Some almost max-length strings for good measure
    kmers.push_back(get_random_dna_string(32, 4));
    kmers.push_back(get_random_dna_string(32, 4));
    kmers.push_back(get_random_dna_string(32, 4));
    kmers.push_back(get_random_dna_string(32, 4));
    kmers.push_back(get_random_dna_string(31, 4));
    kmers.push_back(get_random_dna_string(31, 4));
    kmers.push_back(get_random_dna_string(31, 4));
    kmers.push_back(get_random_dna_string(31, 4));

    // Check all comparisons
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(strings[i] == strings[j]){
                assert(kmers[i] == kmers[j]);
                assert((kmers[i] < kmers[j]) == false);
            } else{
                assert(colex_compare(strings[i], strings[j]) == (kmers[i] < kmers[j]));
            }
        }
    }
}

void test_Kmer(){
    cerr << "Testing k-mer class" << endl;
    Kmer_basic_tests();
    Kmer_colex_tests();
    cerr << "k-mer class tests passed" << endl;
}