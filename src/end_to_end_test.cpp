#include "end_to_end_test.hh"

int main(int argc, char** argv){
    temp_file_manager.set_dir("temp");
    end_to_end_test(argv[1]);
}