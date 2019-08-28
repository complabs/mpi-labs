
CXX=mpic++
CXXFLAGS=-std=c++14 -Wall -O2

ifeq ($(CRAY_PRGENVCRAY), loaded)
   CXX=CC
else ifeq ($(CRAY_PRGENVINTEL), loaded)
   CXX=CC
   CXXFLAGS=$(CXXFLAGS) -D_Float128=__float128
else ifeq ($(CRAY_PRGENVGNU), loaded)
   CXX=CC
endif

all : pi serecv pi-nb psearch golife

clean :
	rm -rf serecv pi pi-nb psearch golife
	rm -rf *.gol core *.save-failed *.out

serecv : serecv.cpp
	$(CXX) $(CXXFLAGS) serecv.cpp -o serecv

pi : pi.cpp
	$(CXX) $(CXXFLAGS) pi.cpp -o pi

pi-nb : pi-nb.cpp
	$(CXX) $(CXXFLAGS) pi-nb.cpp -o pi-nb

psearch : psearch.cpp
	$(CXX) $(CXXFLAGS) psearch.cpp -o psearch

golife : golife.cpp ../lab1/golife.h
	$(CXX) $(CXXFLAGS) golife.cpp -o golife

run : all
	sbatch test-ex.sh 
	qtop
	cat test-ex.out

run-go : all
	mpirun -n 1 ./golife 32 128 200
	mpirun -n 4 ./golife 32 128 200
	paste -d "" life-*.gol
	mpirun -n 8 ./golife 32 128 200
	paste -d "" life-*.gol

verify : golife
	paste -d "" life-*.gol | diff - life.gol
