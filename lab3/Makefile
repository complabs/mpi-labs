
CXX=mpic++
CXXFLAGS=-std=c++14 -Wall -O2
CC=mpicc
CCFLAGS=-Wall

ifeq ($(CRAY_PRGENVCRAY), loaded)
   CXX=CC
   CC=cc
else ifeq ($(CRAY_PRGENVINTEL), loaded)
   CXX=CC
   CC=cc
   CXXFLAGS=$(CXXFLAGS) -D_Float128=__float128
   CCFLAGS=$(CCFLAGS) -D_Float128=__float128
else ifeq ($(CRAY_PRGENVGNU), loaded)
   CXX=CC
   CC=cc
endif

all : dtypes golife-win golife-topology hw-io s1d-topology bw bw-nb latency

clean :
	rm -rf dtypes golife-win golife-topology s1d-topology hw-io bw bw-nb latency
	rm -rf *.gol core *.save-failed out.stl *.out out*.stl

dtypes : dtypes.cpp 
	$(CXX) $(CXXFLAGS) dtypes.cpp -o dtypes

golife-win : golife.cpp ../lab1/golife.h
	$(CXX) $(CXXFLAGS) golife.cpp -o golife-win

golife-topology : golife.cpp ../lab1/golife.h
	$(CXX) $(CXXFLAGS) -DTOPOLOGY golife.cpp -o golife-topology

hw-io : hw-io.cpp
	$(CXX) $(CXXFLAGS) hw-io.cpp -o hw-io

s1d-topology : s1d-topology.cpp
	$(CXX) $(CXXFLAGS) s1d-topology.cpp -o s1d-topology

bw : bw.c
	$(CC) $(CCFLAGS) bw.c -o bw

bw-nb : bw-nb.c
	$(CC) $(CCFLAGS) bw-nb.c -o bw-nb

latency : latency.c
	$(CC) $(CCFLAGS) latency.c -o latency

run : all
	sbatch test-ex.sh 
	qtop
	cat test-ex.out

stl : dtypes
	@echo
	mpirun -n 1 ./dtypes data/sphere.stl out1.stl
	@echo
	@md5sum out1.stl data/sphere.stl
	@echo
	mpirun -n 32 ./dtypes data/sphere.stl out2.stl
	@echo
	@md5sum out2.stl data/sphere.stl
	@echo

