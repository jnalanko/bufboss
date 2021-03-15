.PHONY: tests mergemain buildmain querymain end_to_end_test

all: tests mergemain buildmain querymain end_to_end_test

LIBPATHS=-L sdsl-lite/build/lib/ -L sdsl-lite/build/external/libdivsufsort/lib/ -L stxxl/build/install/lib/
LIBS=-lsdsl -ldivsufsort -ldivsufsort64 -lstxxl -fopenmp -pthread
INCLUDEPATHS=-I sdsl-lite/build/include/ -I sdsl-lite/build/external/libdivsufsort/include -I stxxl/build/install/include/ -I src 

KMC_API_objects=KMC/kmc_api/kmc_file.o KMC/kmc_api/kmer_api.o KMC/kmc_api/mmer.o
KMC_API_includes=KMC/kmc_api

tests:
	mkdir -p bin
	mkdir -p test_out
	$(CXX) -std=c++17 -I src $(LIBPATHS) $(INCLUDEPATHS) src/tests.cpp $(LIBS) -o bin/tests -g -Wall -Wextra -pthread

mergemain:
	$(CXX) -std=c++17 -I src $(LIBPATHS) $(INCLUDEPATHS) src/mergemain.cpp $(LIBS) -o bin/bufboss_update -g -O3 -DNDEBUG -Wall -Wextra -march=native

buildmain:
	$(CXX) -std=c++17 -I src -I $(KMC_API_includes) $(KMC_API_objects) $(LIBPATHS) $(INCLUDEPATHS) src/buildmain.cpp $(LIBS) -o bin/bufboss_build -g -O3 -DNDEBUG -Wall -Wextra -pthread -march=native

querymain:
	$(CXX) -std=c++17 -I src $(LIBPATHS) $(INCLUDEPATHS) src/querymain.cpp $(LIBS) -o bin/bufboss_query -g -O3 -DNDEBUG -Wall -Wextra -march=native

end_to_end_test:
	$(CXX) -std=c++17 -I src $(LIBPATHS) $(INCLUDEPATHS) src/end_to_end_test.cpp $(LIBS) -o bin/end_to_end_test -g -Wall -Wextra -pthread
