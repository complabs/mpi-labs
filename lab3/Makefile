
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

all : golife-win golife-topology hw-io s1d-topology bw bw-nb

clean :
	rm -rf golife-win golife-topology s1d-topology hw-io bw bw-nb
	rm -rf *.gol core *.save-failed

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

run : all
	sbatch test-ex.sh 
	qtop
	cat test-ex.out
