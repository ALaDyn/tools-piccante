/*******************************************************************************
This file is part of tools_pic.

tools_pic is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
tools-pic is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with tools-pic.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <malloc.h>
#include <cmath>
#include <iomanip>
#include <cstring>
#if defined(CINECA)
#include <inttypes.h>
#else
#include <cstdint>
#endif
#include <cstdlib>

void swap_endian(float* in_f, size_t n)
{
  size_t i;
  union {int imio; float fmio; char arr[4];}x;
  char buff;
  for(i=0;i<n;i++)
  {
    x.fmio=in_f[i];
    buff=x.arr[0];
    x.arr[0]=x.arr[3];
    x.arr[3]=buff;
    buff=x.arr[1];
    x.arr[1]=x.arr[2];
    x.arr[2]=buff;
    in_f[i]=x.fmio;
  }
}
void swap_endian(int* in_i, int n)
{
  int i;
  union { int imio; float fmio; char arr[4]; }x;
  char buff;
  for (i = 0; i < n; i++)
  {
    x.imio = in_i[i];
    buff = x.arr[0];
    x.arr[0] = x.arr[3];
    x.arr[3] = buff;
    buff = x.arr[1];
    x.arr[1] = x.arr[2];
    x.arr[2] = buff;
    in_i[i] = x.imio;
  }
}


int is_big_endian();

inline void drawLoadBar(long, long, long, int);

int findIndexMin (double val, float* coords, int numcoords);
int findIndexMax (double val, float* coords, int numcoords);

int main(const int argc, const char *argv[]){
  bool FLAG_lockr[3];
  FLAG_lockr[0] = false;
  FLAG_lockr[1] = false;
  FLAG_lockr[2] = false;

  bool  FLAG_integratex = false,FLAG_integratey = false, FLAG_integratez = false;
  bool FLAG_cutx = false;
  bool FLAG_cuty = false;
  bool FLAG_cutz = false;
  int iminval[] = {0,0,0};
  int imaxval[] = {0,0,0};

  bool FLAG_xmin = false; double xminval = -9999.0;
  bool FLAG_xmax = false; double xmaxval = +9999.0;
  bool FLAG_ymin = false; double yminval = -9999.0;
  bool FLAG_ymax = false; double ymaxval = +9999.0;
  bool FLAG_zmin = false; double zminval = -9999.0;
  bool FLAG_zmax = false; double zmaxval = +9999.0;

  bool doSwap;
  int isFileBigEndian;
  double valueCutx; double valueCuty;
  int Ncells[3], rNproc[3], locOrigin[3], locNcells[3];
  float *xiCoords, *yiCoords, *ziCoords;
  float *riCoords[3];
  int Nproc;
  int Ncomp;
  long size;
  float *fields;
  int* integer_or_halfinteger;
  std::ostringstream nomefile_bin, nomefile_txt;
  nomefile_bin << std::string(argv[1]);
  nomefile_txt << std::string(argv[1]) << ".txt";
  std::ifstream file_bin;
  std::ofstream file_txt;
  file_bin.open(nomefile_bin.str().c_str(), std::ios::binary | std::ios::in);
  file_txt.open(nomefile_txt.str().c_str());
  std::cout << "\nWelcome to the new reader" << std::endl;
  std::cout << "I will read the file: " << nomefile_bin.str() << std::endl;
  if (argc < 1){
    printf("USAGE: reader input_file \n");
  }

  if (file_bin.fail()){
    std::cout << "Input file non trovato" << std::endl;
    return -3;
  }
  for (int i = 2; i < argc; i++){
    if (!std::strncmp(argv[i], "-lockx", 6)){
      FLAG_lockr[0] = true;
    }
    if (!std::strncmp(argv[i], "-locky", 6)){
      FLAG_lockr[1] = true;
    }
    if (!std::strncmp(argv[i], "-lockz", 6)){
      FLAG_lockr[2] = true;
    }
    if (!std::strncmp(argv[i], "-integratex", 11)){
      FLAG_integratex = true;
    }

    if (!std::strncmp(argv[i], "-cutx", 5)){
      valueCutx = atof(argv[i + 1]);
      FLAG_cutx = true;
    }
    if (!std::strncmp(argv[i], "-cuty", 5)){
      valueCuty = atof(argv[i + 1]);
      FLAG_cuty = true;
    }
    if (!std::strncmp(argv[i], "-xmin", 5)){
      FLAG_xmin = true;
      xminval = atof(argv[i + 1]);
    }

    if (!std::strncmp(argv[i], "-xmax", 5)){
      FLAG_xmax = true;
      xmaxval = atof(argv[i + 1]);
    }

    if (!std::strncmp(argv[i], "-ymin", 5)){
      FLAG_ymin = true;
      yminval = atof(argv[i + 1]);
    }

    if (!std::strncmp(argv[i], "-ymax", 5)){
      FLAG_ymax = true;
      ymaxval = atof(argv[i + 1]);
    }

    if (!std::strncmp(argv[i], "-zmin", 5)){
      FLAG_zmin = true;
      zminval = atof(argv[i + 1]);
    }

    if (!std::strncmp(argv[i], "-zmax", 5)){
      FLAG_zmax = true;
      zmaxval = atof(argv[i + 1]);
    }
  }
  if(!(FLAG_lockr[0]||FLAG_lockr[1]||FLAG_lockr[2])){
    FLAG_lockr[2]=true;
  }

  file_bin.read((char*)&isFileBigEndian, sizeof(int));
  doSwap = (isFileBigEndian!=is_big_endian());
  printf("do swap? %i \n", doSwap);
  file_bin.read((char*)Ncells, 3 * sizeof(int));
  if(doSwap)
    swap_endian( Ncells, 3);
  file_bin.read((char*)rNproc, 3 * sizeof(int));
  if(doSwap)
    swap_endian( rNproc, 3);

  Nproc = rNproc[0] * rNproc[1] * rNproc[2];
  xiCoords = new float[Ncells[0]];
  yiCoords = new float[Ncells[1]];
  ziCoords = new float[Ncells[2]];

  riCoords[0] = xiCoords;
  riCoords[1] = yiCoords;
  riCoords[2] = ziCoords;

  file_bin.read((char*)(&Ncomp), sizeof(int));
  if(doSwap)
    swap_endian( &Ncomp, 1);

  for (int c = 0; c < 3; c++){
    file_bin.read((char*)riCoords[c], Ncells[c] * sizeof(float));
    if(doSwap)
      swap_endian( riCoords[c], Ncells[c]);
  }

  std::cout << "IsBigEndian:  " << isFileBigEndian << "\n";
  std::cout << "Ncells:  " << Ncells[0] << "  " << Ncells[1] << "  " << Ncells[2] << "\n";
  std::cout << "Nprocs:  " << rNproc[0] << "  " << rNproc[1] << "  " << rNproc[2] << "\n";
  std::cout << "Ncomp: " << Ncomp << std::endl;
  std::cout << "sizeof long =  " << sizeof(long) << std::endl;

  imaxval[0] = Ncells[0]-1;
  imaxval[1] = Ncells[1]-1;
  imaxval[2] = Ncells[2]-1;

  if (FLAG_xmin)
    iminval[0] = findIndexMin(xminval, riCoords[0], Ncells[0]);

  if (FLAG_xmax)
    imaxval[0] = findIndexMax(xmaxval, riCoords[0], Ncells[0]);

  if (FLAG_ymin)
    iminval[1] = findIndexMin(yminval, riCoords[1], Ncells[1]);

  if (FLAG_ymax)
    imaxval[1] = findIndexMax(ymaxval, riCoords[1], Ncells[1]);

  if (FLAG_zmin)
    iminval[2] = findIndexMin(zminval, riCoords[2], Ncells[2]);

  if (FLAG_zmax)
    imaxval[2] = findIndexMax(zmaxval, riCoords[2], Ncells[2]);

  int allocN[3]={Ncells[0],Ncells[1],Ncells[2]};
  int lockIndex[3]={Ncells[0]/2,Ncells[1]/2,Ncells[2]/2};
  for(int c=0; c<3; c++)
    if(FLAG_lockr[c]){
      allocN[c]=1;
      lockIndex[c]=Ncells[c]/2;

    }
  for(int c=0; c<3; c++)
    printf("FLAG_lockr[%i] = %i  allocN[%i] = %i   lockIndex[%i]=%i \n", c, FLAG_lockr[c], c, allocN[c], c, lockIndex[c]);

  for(int c=0; c < 3; c++)
    size = ((long)Ncomp)*((long)allocN[0])*((long)allocN[1])*((long)allocN[2]);
  fields = new float[size];

  std::cout << "Reading ..." << std::endl; std::cout.flush();


  for (int rank = 0; rank < Nproc; rank++){
    file_bin.read((char*)locOrigin, 3 * sizeof(int));
    if(doSwap)
      swap_endian( locOrigin,3);

    file_bin.read((char*)locNcells, 3 * sizeof(int));
    if(doSwap)
      swap_endian( locNcells,3);

    //        std::cout << "rank= " << rank << "  ";
    //        //std::cout << std::endl;
    //        std::cout << "locNcells: " << locNcells[0] << "  " << locNcells[1] << "  " << locNcells[2] << " ";
    //        std::cout << "orign: " << locOrigin[0] << "  " << locOrigin[1] << "  " << locOrigin[2] << "\n";
    float *locFields;
    int locSize = Ncomp*locNcells[0] * locNcells[1] * locNcells[2];
    locFields = new float[locSize];
    file_bin.read((char*)locFields, locSize*sizeof(float));
    if(doSwap)
      swap_endian( locFields,locSize);

    drawLoadBar(rank + 1, Nproc, Nproc, 30);

    bool flagRead=true;
    for (int c = 0; c < 3; c++){
      if(FLAG_lockr[c])
        if( locOrigin[c]<=lockIndex[c]&& ((locOrigin[c]+locNcells[c])>lockIndex[c]) )
          flagRead = flagRead && true;
        else
          flagRead = false;
    }

    if(flagRead){
      for (int k = 0; k < locNcells[2]; k++){
        for (int j = 0; j < locNcells[1]; j++){
          for (int i = 0; i < locNcells[0]; i++){
            for (int c = 0; c < Ncomp; c++){
              long ii = i + locOrigin[0];
              long jj = j + locOrigin[1];
              long kk = k + locOrigin[2];
              long index = c + Ncomp*ii*(!FLAG_lockr[0]) + Ncomp*allocN[0] * jj* (!FLAG_lockr[1]) + Ncomp*allocN[0] * allocN[1] * kk * (!FLAG_lockr[2]);
              long locIndex = c + Ncomp*i + Ncomp*locNcells[0] * j + Ncomp*locNcells[0] * locNcells[1] * k;
              if(!FLAG_lockr[0] || ii==lockIndex[0])
                if(!FLAG_lockr[1] || jj==lockIndex[1])
                  if(!FLAG_lockr[2] || kk==lockIndex[2])
                    fields[index] = locFields[locIndex];
            }
          }
        }
      }
    }
    delete[] locFields;
  }

  std::cout << std::endl << "Writing to file ..." << std::endl; std::cout.flush();

  {
    if(FLAG_lockr[0]){
      iminval[0]=lockIndex[0];
      imaxval[0]=lockIndex[0]+1;
    }    
    if(FLAG_lockr[1]){     
      iminval[1]=lockIndex[1];
      imaxval[1]=lockIndex[1]+1;
    }    
    if(FLAG_lockr[2]){
      iminval[2]=lockIndex[2];
      imaxval[2]=lockIndex[2]+1;
    }


    for(int c=0; c < 3; c++)
      printf("allocN[%i] = %i   ", c, allocN[c]);
    printf("\n");
    std::stringstream bufstream;

    for (long kk = iminval[2]; kk <imaxval[2] ; kk++){
      for (long jj = iminval[1]; jj <imaxval[1] ; jj++){
        for (long ii = iminval[0]; ii <imaxval[0]; ii++){
          bufstream << std::setw(12) << std::setprecision(5) << xiCoords[ii];
          bufstream << std::setw(12) << std::setprecision(5) << yiCoords[jj];
          bufstream << std::setw(12) << std::setprecision(5) << ziCoords[kk];
          for (int c = 0; c < Ncomp; c++){
            long index = c + Ncomp*ii + Ncomp*allocN[0] * jj + Ncomp*allocN[0] * allocN[1] * kk;
            bufstream << std::setw(12) << std::setprecision(5) << fields[index];
          }
          bufstream << std::endl;
        }
        bufstream << std::endl;
      }
    }
    std::string bufstring = bufstream.str();
    file_txt.write(bufstring.c_str(), bufstring.length());
  }
  std::cout << std::endl;

  file_bin.close();
  file_txt.close();


}


int is_big_endian(){
  union {
    uint32_t i;
    char c[4];
  } bint = { 0x01020304 };

  return bint.c[0] == 1;
}

inline void drawLoadBar(long i, long Ntot, long R, int sizeBar){
  if (i % (Ntot / R) != 0) return;

  std::cout << "\r";

  float ratio = i / ((float)Ntot);

  int numSymbols = (int)(sizeBar*ratio);

  std::cout << std::setw(3) << (int)(ratio * 100) << "% [";

  int j;
  for (j = 0; j < numSymbols; j++){
    std::cout << "=";
  }
  for (; j < sizeBar; j++){
    std::cout << " ";
  }

  std::cout << "]";
  std::cout.flush();
}


int findIndexMin (double val, float* coords, int numcoords){
  if (numcoords <= 1)
    return 0;

  if (val <= coords[0])
    return 0;

  for (int i = 1; i < numcoords; i++){
    if (val < coords[i])
      return (i-1);
  }

}

int findIndexMax (double val, float* coords, int numcoords){
  if (numcoords<= 1)
    return 0;

  if (val >= coords[numcoords-1])
    return (numcoords-1);

  for (int i = (numcoords-1); i >= 0; i--){
    if (val > coords[i])
      return (i+1);
  }


}

