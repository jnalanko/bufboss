#include "BOSS/BOSS_tests.hh"
#include "test_Kmer.hh"
#include "test_light_merging.hh"
#include "end_to_end_test.hh"
#include "test_boss_builder.hh"

int main(){
    srand(1234);
    temp_file_manager.set_dir("./temp");
    test_BOSS_builder();
    test_light_merging<BOSS<sdsl::bit_vector>>();
    test_Kmer();
    test_BOSS<BOSS<sdsl::bit_vector>>();
    //write_log("Warning: some tests disabled");
    end_to_end_test("data/0.fna");
}
