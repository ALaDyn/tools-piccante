COMPILER = g++
EXE = piccante

OPT = -O3 -std=c++0x

LFLAGS = -Wall

LIB = 
EXE1 = new_reader
EXE2 = WangularSpectrum
EXE3 = WenergySpectrum
EXE4 = WsliceDaBinario
EXE5 = light_reader

SRC1 = newfield_reader.cpp
SRC2 = WangularSpectrum.cpp
SRC3 = WenergySpectrum.cpp
SRC4 = WsliceDaBinario.cpp
SRC5 = light-reader.cpp

all : $(EXE1) $(EXE2) $(EXE3) $(EXE4)  $(EXE5) 

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

clean :
				rm $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5)


