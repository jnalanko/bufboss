#pragma once

#include "Light_merging.hh"
#include "BOSS/BOSS_tests.hh"
#include "BufBOSS.hh"

template<typename boss_t = BOSS<sdsl::bit_vector>>
class DynBufferMergeLightTester{
public:
    void test_light_merging(){
        Boss_Tester tester;
        auto testcases = tester.generate_testcases();;
        for(LL testcase_id = 0; testcase_id < (LL)testcases.size(); testcase_id++){
            set<string> true_edgemers;
            if(testcases[testcase_id].k+1 > 32) continue; // Not supported by the packed k-mer representation
            LL k = testcases[testcase_id].k;
            BufBOSS<boss_t> bufboss(k);
            vector<string> reads;
            for(string read : testcases[testcase_id].reads){
                reads.push_back(read);

                bufboss.add_edgemers_to_buffer_not_including_dummies(reads.back());
                bufboss.flush();
                
                for(string edgemer : get_all_distinct_kmers(read, k+1))
                    true_edgemers.insert(edgemer);

                set<string> our_edgemers = bufboss.boss.get_full_edgemers();
                assert(our_edgemers == true_edgemers);
            }
        }
    }

    void colex_print(const set<string>& S){
        vector<string> v(S.begin(), S.end());
        sort(v.begin(), v.end(), colex_compare);
        cout << v << endl;
    }

    // Adds the read to bufboss and to truth
    void add_read(BufBOSS<boss_t>& bufboss, set<string>& true_edgemers, string read){
        bufboss.add_edgemers_to_buffer_not_including_dummies(read);
        for(string edgemer : get_all_distinct_kmers(read, bufboss.boss.get_k()+1)) 
            true_edgemers.insert(edgemer);
    }

    // Deletes the read from bufboss and from truth
    void del_read(BufBOSS<boss_t>& bufboss, set<string>& true_edgemers, string read){
        bufboss.delete_edgemers_from_buffer_not_including_dummies(read);
        for(string edgemer : get_all_distinct_kmers(read, bufboss.boss.get_k()+1)) 
            if(true_edgemers.count(edgemer)) true_edgemers.erase(edgemer);
    }

    void test_interleaved_additions_and_deletions(){
        Boss_Tester tester;
        auto testcases = tester.generate_testcases();;
        for(LL testcase_id = 0; testcase_id < (LL)testcases.size(); testcase_id++){
            if(testcases[testcase_id].k+1 > 32) continue; // Not supported by the packed k-mer representation
            LL k = testcases[testcase_id].k;
            BufBOSS<boss_t> bufboss(k);
            vector<string> reads = testcases[testcase_id].reads;
            set<string> true_edgemers;
            for(LL round = 0; round < (LL)reads.size(); round++){
                cout << tester.BOSS_to_string(bufboss.boss) << endl;
                cout << true_edgemers << endl;
                del_read(bufboss, true_edgemers, reads[rand()%reads.size()]);
                cout << true_edgemers << endl;
                add_read(bufboss, true_edgemers, reads[rand()%reads.size()]);
                cout << true_edgemers << endl;
                del_read(bufboss, true_edgemers, reads[rand()%reads.size()]);
                cout << true_edgemers << endl;
                add_read(bufboss, true_edgemers, reads[rand()%reads.size()]);
                cout << true_edgemers << endl;

                bufboss.flush();

                // Check that the set of (k+1)-mers is what is should be
                set<string> our_edgemers = bufboss.boss.get_full_edgemers();
                assert(our_edgemers == true_edgemers);
                assert(bufboss.boss.has_exactly_one_source_node());
            }
        }
    }

    void test_permute_boss_labels(){
        cerr << "Permuting BOSS labels test" << endl;

        Boss_Tester<boss_t> tester;
        
        for(typename Boss_Tester<boss_t>::TestCase tcase : tester.generate_testcases()){
            
            if(tcase.k+1 > 32) continue; // Not supported by the packed k-mer representation
            cout << "testcase" << endl;
            LL k = tcase.k;
            boss_t boss = build_BOSS_with_maps(tcase.reads, k);
            BufBOSS<boss_t> bufboss(boss);
            DynBufferMergeLight<boss_t> merger(bufboss.boss, bufboss.addbuf, bufboss.delbuf);
            //string outlabels = boss.get_outlabels();
            Packed_string outlabels({'A','C','G','T'}, boss.number_of_edges(), 'A');
            for(LL i = 0; i < boss.number_of_edges(); i++)
                outlabels[i] = boss.outlabels_at(i);
            
            Packed_string labels_in({'\0','A','C','G','T'}, boss.number_of_nodes(), '\0');
            for(LL node = 0; node < boss.number_of_nodes(); node++){
                labels_in[node] = boss.incoming_character(node);
            }

            Packed_string labels_out({'\0','A','C','G','T'}, labels_in.size(), '\0');
            vector<string> rev_kmers(labels_in.size());
            for(LL i = 0; i < k; i++){
                for(LL j = 0; j < (LL)labels_in.size(); j++){
                    if(labels_in[j] != '\0') rev_kmers[j] += labels_in[j];
                }
                merger.permute_boss_labels(labels_in, labels_out, outlabels);
                labels_in = labels_out;
            }
            for(LL node = 0; node < boss.number_of_nodes(); node++){
                char truth[k+1];
                boss.get_node_label(node, truth);
                string kmer(rev_kmers[node].rbegin(), rev_kmers[node].rend());
                assert(strcmp(truth, kmer.data()) == 0);
            }
        }

        cerr << "Permuting BOSS labels test passed" << endl;
    }

};

void test_Packed_string(){
    LL n = 10;
    Packed_string vec({'A','C','G','T','X'}, n, 'C');
    assert(vec.size() == n);

    assert(vec.log2_ceiling(1) == 0);
    assert(vec.log2_ceiling(2) == 1);
    assert(vec.log2_ceiling(3) == 2);
    assert(vec.log2_ceiling(4) == 2);
    assert(vec.log2_ceiling(5) == 3);
    assert(vec.log2_ceiling(6) == 3);
    assert(vec.log2_ceiling(7) == 3);
    assert(vec.log2_ceiling(8) == 3);

    for(LL i = 0; i < n; i++) assert(vec[i] == 'C');
    vec[0] = 'A';
    vec[1] = 'C';
    vec[2] = 'G';
    vec[3] = 'T';
    vec[4] = 'X';
    assert(vec[0] == 'A');
    assert(vec[1] == 'C');
    assert(vec[2] == 'G');
    assert(vec[3] == 'T');
    assert(vec[4] == 'X');
    for(LL i = 0; i < n; i++) cout << vec[i];
    cout << endl;

    Packed_string vec2({'A','C','G','T'}, n,'T');
    vec[4] = vec2[6]; // # vec2[6] returns a char_ref, which is cast to char implicitly and assigned to the reference at vec[4]
    for(LL i = 0; i < n; i++) cout << vec[i];
    cout << endl;
    assert(vec[4] == vec2[6]);
    
    cerr << "Packed_string test passed" << endl;
}

template<typename boss_t>
void test_light_merging(){
    test_Packed_string();
    DynBufferMergeLightTester<boss_t> tester;
    tester.test_interleaved_additions_and_deletions();
    tester.test_permute_boss_labels();
    tester.test_light_merging();
}