DTL_BASE := ./bloomfilter-bsd
DTL_FLAGS := -Wno-ignored-attributes -Wno-unused-variable -Wno-sign-compare -DUSE_DTL
DTL_LIBFILE := dtl_core-avx2
DTL_CC := $(shell which gcc)
DTL_CXX := $(shell which g++)


# WHEN BUILDING FOR LINKING INTO ROCKSDB EXPERIMENTS, DO DECLARE IN ROCKSDB MAKEFILE TOO
DTL_SUFFIX := -DKEEP_STATS

ifeq ($(use_dtl),1)
	DTL_SUFFIX := -I$(DTL_BASE)/src -I$(DTL_BASE)/thirdparty -L$(DTL_BASE) -l$(DTL_LIBFILE) $(DTL_FLAGS) -lpthread -lnuma
endif

DSTOPT=-Ofast
#DSTOPT=-g

all: test # run_workload bloom

#surf.so: ../SuRF_standalone/SuRF/include/surf.hpp
#	g++ -Ofast -march=native -fpic -shared -x c++ ../SuRF_standalone/SuRF/include/surf.hpp -o surf.so

#dst.h.gch: dst.h $(DTL_BASE)/lib$(DTL_LIBFILE).a
#	g++ -Ofast -march=native -c dst.h -o dst.h.gch $(DTL_SUFFIX) -Wall -Wextra #-I../SuRF_standalone/SuRF/include/

dst.o: dst.cpp dst.h MurmurHash3.o util.h
	g++ $(DSTOPT) -march=native -c dst.cpp --std=c++17 $(DTL_SUFFIX) -Wall -Wextra -I../SuRF_standalone/SuRF/include -g  -Wl,--stack=512000000

MurmurHash3.o: MurmurHash3.cpp MurmurHash3.h
	g++ -Ofast -march=native -c MurmurHash3.cpp

test: MurmurHash3.o dst.o test.cpp
	g++ -Ofast -march=native MurmurHash3.o dst.o test.cpp -o test --std=c++17 $(DTL_SUFFIX)

string_test: MurmurHash3.o dst.o string_test.cpp
	g++ -Ofast -march=native MurmurHash3.o dst.o string_test.cpp -o string_test --std=c++17
	
bloom: MurmurHash3.o dst.o bloom.cpp 
	g++ -Ofast -march=native MurmurHash3.o dst.o bloom.cpp -o bloom --std=c++17 $(DTL_SUFFIX)

run_workload: MurmurHash3.o dst.o run_workload.cpp
	g++ -Ofast -march=native MurmurHash3.o dst.o run_workload.cpp -o run_workload --std=c++17 $(DTL_SUFFIX) -Wall -Wextra -g

run_workload_string: MurmurHash3.o dst.o run_workload.cpp
	g++ -Ofast -march=native MurmurHash3.o dst.o run_workload.cpp -o run_workload_string --std=c++17 $(DTL_SUFFIX) -Wall -Wextra -lpthread -DKEYTYPE=string -g

run_workload_disk: MurmurHash3.o dst.o run_workload.cpp disk.h
	g++ -Ofast -march=native MurmurHash3.o dst.o run_workload.cpp -o run_workload_disk --std=c++17 $(DTL_SUFFIX) -Wall -Wextra -lpthread -DUSE_DISK -g

run_workload_string_disk: MurmurHash3.o dst.o run_workload.cpp disk.h
	g++ -Ofast -march=native MurmurHash3.o dst.o run_workload.cpp -o run_workload_string_disk --std=c++17 $(DTL_SUFFIX) -Wall -Wextra -lpthread -DUSE_DISK -DKEYTYPE=string -g

disk_experiment: MurmurHash3.o dst.o disk_experiment.cpp
	g++ -Ofast -march=native MurmurHash3.o dst.o disk_experiment.cpp -o disk_experiment --std=c++17 $(DTL_SUFFIX) -Wall -Wextra


$(DTL_BASE)/lib$(DTL_LIBFILE).a:
	cd $(DTL_BASE) && export CC=$(DTL_CC) export CXX=$(DTL_CXX) && cmake . && make $(DTL_LIBFILE)

dtltest: dtltest.cpp $(DTL_BASE)/lib$(DTL_LIBFILE).a
	g++ -O0 -g -march=native MurmurHash3.o dst.o dtltest.cpp -o dtltest --std=c++17 $(DTL_SUFFIX)

clean:
	rm -f *.o *.h.gch test run_workload bloom dtltest
