
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


all : hw serecv pi psearch golife

clean :
	rm -rf hw serecv pi psearch golife
	rm -rf *.gol core *.save-failed *.out

hw : hw.cpp
	$(CXX) $(CXXFLAGS) hw.cpp -o hw

serecv : serecv.cpp
	$(CXX) $(CXXFLAGS) serecv.cpp -o serecv

pi : pi.cpp
	$(CXX) $(CXXFLAGS) pi.cpp -o pi

psearch : psearch.cpp
	$(CXX) $(CXXFLAGS) psearch.cpp -o psearch

golife : golife.cpp golife.h
	$(CXX) $(CXXFLAGS) golife.cpp -o golife

run : all
	sbatch test-ex.sh 
	qtop
	cat test-ex.out

run-go : all
	mpirun -n 4 ./golife 100 100 2000
	@echo "-----------------------------------------------------------------------------"
	@cat life-0.gol
	@echo "-----------------------------------------------------------------------------"
	@cat life-1.gol
	@echo "-----------------------------------------------------------------------------"
	@cat life-2.gol
	@echo "-----------------------------------------------------------------------------"
	@cat life-3.gol
	@echo "-----------------------------------------------------------------------------"

verify : golife
	cat life-*.gol | diff - life.gol
