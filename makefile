MPICOMPILER = mpic++
OPT = -O3 -std=c++11
MPI_OPT = -O3 -std=c++11
SRC_FOLDER = src

FILES = newReader.cpp \
WangularSpectrum.cpp \
WenergySpectrum.cpp \
WsliceDaBinario.cpp \
WspacialSpectrum.cpp \
lightReader.cpp \
frogReader.cpp \
titan.cpp \
titan2.cpp \
clusters.cpp \
multiFrogReader.cpp

MPIFILES = MPItitan.cpp \
MPIpicker.cpp

DEPFILES = utilities-tools.cpp

SRC = $(addprefix $(SRC_FOLDER)/, $(FILES))
MPI_SRC = $(addprefix $(SRC_FOLDER)/, $(MPIFILES))
DEP_SRC = $(addprefix $(SRC_FOLDER)/, $(DEPFILES))

MPIEXES = $(addsuffix .mpi, $(basename $(MPIFILES)))
EXES = $(addsuffix .exe, $(basename $(FILES)))
DEPS = $(addsuffix .o, $(basename $(DEPFILES)))

all : $(DEPS)
all : $(EXES)
all : $(MPIEXES)

boost : OPT = -O3 -DUSE_BOOST
boost : LIB = -lboost_random 
boost : $(DEPS)
boost : $(EXES)
boost : $(MPIEXES)

fermi : OPT = -O3 -std=c++0x
fermi : MPICOMPILER = mpixlcxx
fermi : MPI_OPT = -qipa=partition=large -qarch=qp -qtune=qp -qmaxmem=-1
fermi : $(DEPS)
fermi : $(EXES)
fermi : $(MPIEXES)

cnaf : OPT = -O3 -std=c++0x 
cnaf : MPI_OPT = -O3 -std=c++0x 
cnaf : $(DEPS)
cnaf : $(EXES)
cnaf : $(MPIEXES)

debug : OPT = -O0 -g -std=c++0x 
debug : MPI_OPT = -O0 -g -std=c++0x 
debug : $(DEPS)
debug : $(EXES)
debug : $(MPIEXES)

%.o: $(SRC_FOLDER)/%.cpp
	$(CXX) $(OPT) -c -o $@ $<

%.exe: $(SRC_FOLDER)/%.cpp
	$(CXX) $(OPT) -o $@ $(DEPS) $< $(LIB)

%.mpi: $(SRC_FOLDER)/%.cpp
	$(MPICOMPILER) $(MPI_OPT) -o $@ $< $(LIB)

clean:
	rm -f $(DEPS)

cleanall: clean
	rm -f $(EXES) $(MPIEXES)


