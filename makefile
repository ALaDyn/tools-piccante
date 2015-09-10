COMPILER = g++
MPICOMPILER = mpic++
OBJ = utilities-tools.o

OPT = -O3 -std=c++0x


LFLAGS = -Wall

LIB = 
EXE1 = newReader
EXE2 = WangularSpectrum
EXE3 = WenergySpectrum
EXE4 = WsliceDaBinario
EXE5 = lightReader
EXE6 = WspacialSpectrum
EXE7 = titan
EXE8 = MPItitan
Exe9 = frogReader
EXE10 = clusters
EXE11 = titan2
EXE12 = multiFrogReader
EXE13 = MPIpicker


SRC1 = newReader.cpp
SRC2 = WangularSpectrum.cpp
SRC3 = WenergySpectrum.cpp
SRC4 = WsliceDaBinario.cpp
SRC5 = lightReader.cpp
SRC6 = WspacialSpectrum.cpp
SRC7 = titan.cpp
SRC8 = MPItitan.cpp
SRC9 = frogReader.cpp
SRC10 = clusters.cpp
SRC11 = titan2.cpp
SRC12 = multiFrogReader.cpp
SRC13 = MPIpicker.cpp


EXES = $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10)  $(EXE11)  $(EXE12) $(EXE13)
all : $(EXES)

boost : OPT = -O3 -DUSE_BOOST
boost : LIB = -lboost_random 
boost: $(EXES)
debug : OPT = -O0 -g 
debug : $(EXE)


$(EXE1) : $(SRC1)
	$(COMPILER) $(SRC1) -o $(EXE1) $(OPT)   $(LIB) 

$(EXE2) : $(SRC2)
	$(COMPILER) $(SRC2) -o $(EXE2) $(OPT)   $(LIB) 

$(EXE3) : $(SRC3)
	$(COMPILER) $(SRC3) -o $(EXE3) $(OPT)   $(LIB) 

$(EXE4) : $(SRC4)
	$(COMPILER) $(SRC4) -o $(EXE4) $(OPT)   $(LIB) 

$(EXE5) : $(SRC5)
	$(COMPILER) $(SRC5) -o $(EXE5) $(OPT)   $(LIB)

$(EXE6) : $(SRC6)
	$(COMPILER) $(SRC6) -o $(EXE6) $(OPT)   $(LIB)

$(EXE7) : $(SRC7)
	$(COMPILER) $(SRC7) -o $(EXE7) $(OPT)   $(LIB)

$(EXE8) : $(SRC8)
	$(MPICOMPILER) $(SRC8) -o $(EXE8) $(OPT) $(LIB)

$(EXE9) : $(SRC9)
	$(COMPILER) $(SRC9) -o $(EXE9) $(OPT)   $(LIB)

utilities-tools.o: utilities-tools.cpp
	$(COMPILER)  $(OPT) -c utilities-tools.cpp
$(EXE10) : $(SRC10) utilities-tools.o
	$(COMPILER) $(SRC10) -o $(EXE10) $(OPT) $(OBJ) $(LIB)
$(EXE11) : $(SRC11)
	$(COMPILER) $(SRC11) -o $(EXE11) $(OPT)   $(LIB)

$(EXE12) : $(SRC12)
	$(COMPILER) $(SRC12) -o $(EXE12) $(OPT)   $(LIB)

$(EXE13) : $(SRC13)
	$(MPICOMPILER) $(SRC13) -o $(EXE13) $(OPT) $(LIB)

clean :
	rm $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(EXE13)  utilities-tools.o


