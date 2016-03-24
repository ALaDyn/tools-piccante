
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
#include "utilities-tools.h"


struct ALL_FLAGS{
    bool FLAG_xmin, FLAG_ymin, FLAG_zmin;
    bool FLAG_xmax, FLAG_ymax, FLAG_zmax;
    double xminval, yminval, zminval;
    double xmaxval, ymaxval, zmaxval;
    double lockValue[3];
    bool FLAG_lockr[3];
    int iminval[3];
    int imaxval[3];
    int sampling[3];
    bool vtkOutput;
    bool doSwap;

    ALL_FLAGS(){
      FLAG_xmin = FLAG_ymin = FLAG_zmin = false;
      FLAG_xmax = FLAG_ymax = FLAG_zmax = false;
      xminval = yminval = zminval = -9999.0;
      xmaxval = ymaxval = zmaxval = +9999.0;
      FLAG_lockr[0] = false;
      FLAG_lockr[1] = false;
      FLAG_lockr[2] = false;
      iminval[0] = iminval[1] = iminval[2] = 0;
      imaxval[0] = imaxval[1] = imaxval[2] = 0;
      sampling[0] = sampling[1] = sampling[2] = 1;
      vtkOutput = false;

    }
};
struct OUTPUT_DATA{
    int lockIndex[3];
    int allocN[3];
    uint64_t size;
    uint64_t NNN;
    float *savedFields;
};

struct FILE_DATA{
    int isBigEndian;
    int dataNSize[3], rNproc[3];

    float *riCoords[3];
    int Nproc;
    int Ncomp;

};

std::string composeFileName(std::string strippedFileName, int fileId);
void message(std::string msg);
void errorMessage(std::string msg);
int howManyFilesExist(std::string strippedFileName);
int getAndCheckFileNumber(std::string strippedFileName);

inline void drawLoadBar(uint64_t, uint64_t, uint64_t, int);

int findIndexMin (double val, float* coords, int numcoords);
int findIndexMax (double val, float* coords, int numcoords);
void parseInputArguments(ALL_FLAGS &myFlags, int argc, char **argv);
void writeUsageHelper();
void checkFlags(ALL_FLAGS &myFlags, FILE_DATA &fileData, OUTPUT_DATA &outputData);

void updateFileIfEOF(std::ifstream &file_bin, int &fileId, std::string strippedName);

void readFromFile(std::ifstream &file_bin, int *buffer, int N, bool doSwap);
void readFromFile(std::ifstream &file_bin, float *buffer, int N, bool doSwap);


void printVTKFile(OUTPUT_DATA &outputData, ALL_FLAGS &myFlags, FILE_DATA fileData, std::string outputfileName);


int main( int narg,  char **args){
  ALL_FLAGS myFlags;
  FILE_DATA fileData;
  OUTPUT_DATA outputData;
  bool doSwap;


  std::cout<< "  There  are " << narg << " arguments" << std::endl;

  if (narg < 2){
    writeUsageHelper();
    exit(0);
  }


  std::ostringstream inputString, outputfileName;
  inputString << std::string(args[1]);
  std::ifstream file_bin;
  int fileNumber=getAndCheckFileNumber(inputString.str());
  std::cout << "\nWelcome to the new reader" << std::endl;
  std::cout<< "  There seems to be " << fileNumber << " files" << std::endl;



  parseInputArguments(myFlags, narg, args);

  int fileId = 0;
  std::cout << "the input string is >>> " << inputString.str() << std::endl;
  std::string currentFileName = composeFileName(inputString.str(),fileId);
  std::cout << "... opening the file >>> " << currentFileName << std::endl;
  file_bin.open(currentFileName.c_str(), std::ios::binary | std::ios::in);
  if (file_bin.fail()){
    std::cout << "Input file non trovato" << std::endl;
    return -3;
  }



  file_bin.read((char*)&fileData.isBigEndian, sizeof(int));
  doSwap = (fileData.isBigEndian!=is_big_endian());
  printf("do swap? %i \n", doSwap);

  file_bin.read((char*)fileData.dataNSize, 3 * sizeof(int));
  if(doSwap)
    swap_endian( fileData.dataNSize, 3);
  file_bin.read((char*)fileData.rNproc, 3 * sizeof(int));
  if(doSwap)
    swap_endian( fileData.rNproc, 3);

  fileData.Nproc = fileData.rNproc[0] * fileData.rNproc[1] * fileData.rNproc[2];
  fileData.riCoords[0] = new float[fileData.dataNSize[0]];
  fileData.riCoords[1] = new float[fileData.dataNSize[1]];
  fileData.riCoords[2] = new float[fileData.dataNSize[2]];

  file_bin.read((char*)(&fileData.Ncomp), sizeof(int));
  if(doSwap)
    swap_endian( &fileData.Ncomp, 1);

  for (int c = 0; c < 3; c++){
    file_bin.read((char*)fileData.riCoords[c], fileData.dataNSize[c] * sizeof(float));
    if(doSwap)
      swap_endian( fileData.riCoords[c], fileData.dataNSize[c]);
  }

  std::cout << "fileData.isBigEndian:  " << fileData.isBigEndian << "\n";
  std::cout << "Ncells:  " << fileData.dataNSize[0] << "  " << fileData.dataNSize[1] << "  " << fileData.dataNSize[2] << "\n";
  std::cout << "Nprocs:  " << fileData.rNproc[0] << "  " << fileData.rNproc[1] << "  " << fileData.rNproc[2] << "\n";
  std::cout << "Ncomp: " << fileData.Ncomp << std::endl;
  std::cout << "sizeof uint64_t =  " << sizeof(uint64_t) << std::endl;



  checkFlags(myFlags, fileData, outputData);


  outputData.size = ((uint64_t)fileData.Ncomp)*((uint64_t)outputData.allocN[0])*((uint64_t)outputData.allocN[1])*((uint64_t)outputData.allocN[2]);
  outputData.savedFields = new float[outputData.size];




  for(int c=0; c<3; c++)
    printf("FLAG_lockr[%i] = %i  outputData.allocN[%i] = %i   outputData.lockIndex[%i]=%i \n", c, myFlags.FLAG_lockr[c], c, outputData.allocN[c], c, outputData.lockIndex[c]);


  std::cout << "Reading ..." << std::endl; std::cout.flush();


  for (int rank = 0; rank < fileData.Nproc; rank++){

    updateFileIfEOF(file_bin, fileId, inputString.str());

    int locOrigin[3], locNcells[3];

    readFromFile(file_bin, locOrigin, 3, doSwap);
    readFromFile(file_bin, locNcells, 3, doSwap);

    float *localFields;
    int locSize = fileData.Ncomp*locNcells[0] * locNcells[1] * locNcells[2];
    localFields = new float[locSize];

    readFromFile(file_bin, localFields, locSize, doSwap);


    drawLoadBar(rank + 1, fileData.Nproc, fileData.Nproc, 30);

    bool flagRead=true;
    for (int c = 0; c < 3; c++){
      if(myFlags.FLAG_lockr[c])
        if( locOrigin[c]<=outputData.lockIndex[c]&& ((locOrigin[c]+locNcells[c])>outputData.lockIndex[c]) )
          flagRead = flagRead && true;
        else
          flagRead = false;
    }
    int shouldWrite[3];
    if(flagRead){
      int savedI[3], globalI[3];
      for (int k = 0; k < locNcells[2]; k++){
        globalI[2] = k + locOrigin[2];
        savedI[2] = (globalI[2] - myFlags.iminval[2])/myFlags.sampling[2];
        shouldWrite[2] = !((globalI[2] - myFlags.iminval[2])%myFlags.sampling[2]);
        for (int j = 0; j < locNcells[1]; j++){
          globalI[1] = j + locOrigin[1];
          savedI[1] = (globalI[1] - myFlags.iminval[1])/myFlags.sampling[1];
          shouldWrite[1] = !((globalI[1] - myFlags.iminval[1])%myFlags.sampling[1]);

          for (int i = 0; i < locNcells[0]; i++){
            globalI[0] = i + locOrigin[0];
            savedI[0] = (globalI[0] - myFlags.iminval[0])/myFlags.sampling[0];
            shouldWrite[0] = !((globalI[0] - myFlags.iminval[0])%myFlags.sampling[0]);

            for (uint64_t c = 0; c < fileData.Ncomp; c++){
              uint64_t index = c*outputData.allocN[0]*outputData.allocN[1]*outputData.allocN[2] + savedI[0] + outputData.allocN[0] * savedI[1] +  outputData.allocN[0]*outputData.allocN[1]*savedI[2];
              uint64_t locIndex = c + fileData.Ncomp*i + fileData.Ncomp*locNcells[0] * j + fileData.Ncomp*locNcells[0] * locNcells[1] * k;

              if(globalI[0] >= myFlags.iminval[0] && globalI[0] < myFlags.imaxval[0] && shouldWrite[0])
                if(globalI[1] >= myFlags.iminval[1] && globalI[1] < myFlags.imaxval[1] && shouldWrite[1])
                  if(globalI[2] >= myFlags.iminval[2] && globalI[2] < myFlags.imaxval[2] && shouldWrite[2]){
                    if (index >= outputData.size){
                      std::cout<< index << "( " << outputData.size << " ) " << c << " " << i << " " << j << " " << k << " " << savedI[0] << " " << savedI[1] << " " << savedI[2] << " " << std::endl;
                    }
                    else
                      outputData.savedFields[index] = localFields[locIndex];
                  }
            }
          }
        }
      }
    }
    delete[] localFields;
  }
  file_bin.close();
  std::cout << std::endl << "Writing to file ..." << std::endl; std::cout.flush();

  if(myFlags.vtkOutput){

    printf("VTK enabled\n");
    outputfileName << std::string(args[1]) << ".vtk";
    printVTKFile(outputData, myFlags, fileData, outputfileName.str());


  }
  else{
    outputfileName << std::string(args[1]) << ".txt";
    std::ofstream file_txt;
    uint64_t totPts = outputData.allocN[0] * outputData.allocN[1]* outputData.allocN[2];
    file_txt.open(outputfileName.str().c_str());
    for(int c=0; c < 3; c++)
      printf("allocN[%i] = %i   ", c, outputData.allocN[c]);
    printf("\n");
    std::stringstream bufstream;
    for (uint64_t kk = 0; kk <outputData.allocN[2] ; kk++){
      for (uint64_t jj = 0; jj <outputData.allocN[1] ; jj++){
        for (uint64_t ii = 0; ii <outputData.allocN[0]; ii++){
          int global[3];
          global[0] = ii*myFlags.sampling[0] + myFlags.iminval[0];
          global[1] = jj*myFlags.sampling[1] + myFlags.iminval[1];
          global[2] = kk*myFlags.sampling[2] + myFlags.iminval[2];

          drawLoadBar(ii + (jj)*outputData.allocN[0] + kk*outputData.allocN[0]*outputData.allocN[1] + 1, totPts, outputData.allocN[0], 30);

          bufstream << std::setw(12) << std::setprecision(5) << fileData.riCoords[0][global[0]];
          bufstream << std::setw(12) << std::setprecision(5) << fileData.riCoords[1][global[1]];
          bufstream << std::setw(12) << std::setprecision(5) << fileData.riCoords[2][global[2]];
          for (int c = 0; c < fileData.Ncomp; c++){
            //uint64_t index = c + fileData.Ncomp*ii + fileData.Ncomp*outputData.allocN[0] * jj + fileData.Ncomp*outputData.allocN[0] * outputData.allocN[1] * kk;
            uint64_t index = c*outputData.allocN[0]*outputData.allocN[1]*outputData.allocN[2] + ii + outputData.allocN[0]*jj +  outputData.allocN[0]*outputData.allocN[1]*kk;

            bufstream << std::setw(12) << std::setprecision(5) << outputData.savedFields[index];
          }
          bufstream << std::endl;
        }
        bufstream << std::endl;
      }
    }
    std::string bufstring = bufstream.str();
    file_txt.write(bufstring.c_str(), bufstring.length());
    file_txt.close();
  }
  std::cout << std::endl;




}



inline void drawLoadBar(uint64_t i, uint64_t Ntot, uint64_t R, int sizeBar){
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

std::string composeFileName(std::string strippedFileName, int fileId){
  const int numZeroes = 3;
  std::stringstream ss;
  ss << strippedFileName <<"."<< std::setfill('0') << std::setw(numZeroes) << fileId;

  return ss.str();
}

void message(std::string msg){
  std::cout << msg << std::endl;
  std::cout.flush();
  return;
}

void errorMessage(std::string msg){
  std::cout << "ERROR: "  << msg << std::endl;
  std::cout.flush();
  exit(1);

}


int howManyFilesExist(std::string strippedFileName){
  int counter = 0;
  int fID = 0;
  bool dummyFlag = true;

  while(dummyFlag){
    std::string bstring = composeFileName(strippedFileName, fID);
    std::ifstream infile;
    infile.open(bstring.c_str());
    if(infile.good()){
      counter++;
      fID++;
      infile.close();
    }
    else{
      dummyFlag=false;
    }
  }
  if(counter == 0){
    std::ifstream infile;
    infile.open(strippedFileName.c_str());
    if(infile.good()){
      counter++;
      infile.close();
      //single_file_wo_suffix = true;
    }
  }
  return counter;
}



int getAndCheckFileNumber(std::string strippedFileName){
  int fileNumber = howManyFilesExist(strippedFileName);
  if(fileNumber <= 0)
    errorMessage("Input file not found");

  return fileNumber;
}

void parseInputArguments(ALL_FLAGS &myFlags, int argc, char **argv){
  for (int i = 2; i < argc; i++){
    if (!std::strncmp(argv[i], "-lockx", 6)){
      myFlags.FLAG_lockr[0] = true;
      myFlags.lockValue[0]= atof(argv[i + 1]);
    }
    if (!std::strncmp(argv[i], "-locky", 6)){
      myFlags.FLAG_lockr[1] = true;
      myFlags.lockValue[1]= atof(argv[i + 1]);
    }
    if (!std::strncmp(argv[i], "-lockz", 6)){
      myFlags.FLAG_lockr[2] = true;
      myFlags.lockValue[2]= atof(argv[i + 1]);
    }

    if (!std::strncmp(argv[i], "-xsample", 8)){
      myFlags.sampling[0] = atoi(argv[i + 1]);;
    }
    if (!std::strncmp(argv[i], "-ysample", 8)){
      myFlags.sampling[1] = atoi(argv[i + 1]);;
    }
    if (!std::strncmp(argv[i], "-zsample", 8)){
      myFlags.sampling[2] = atoi(argv[i + 1]);;
    }

    if (!std::strncmp(argv[i], "-xmin", 5)){
      myFlags.FLAG_xmin = true;
      myFlags.xminval = atof(argv[i + 1]);
    }

    if (!std::strncmp(argv[i], "-xmax", 5)){
      myFlags.FLAG_xmax = true;
      myFlags.xmaxval = atof(argv[i + 1]);
    }

    if (!std::strncmp(argv[i], "-ymin", 5)){
      myFlags.FLAG_ymin = true;
      myFlags.yminval = atof(argv[i + 1]);
    }

    if (!std::strncmp(argv[i], "-ymax", 5)){
      myFlags.FLAG_ymax = true;
      myFlags.ymaxval = atof(argv[i + 1]);
    }

    if (!std::strncmp(argv[i], "-zmin", 5)){
      myFlags.FLAG_zmin = true;
      myFlags.zminval = atof(argv[i + 1]);
    }

    if (!std::strncmp(argv[i], "-zmax", 5)){
      myFlags.FLAG_zmax = true;
      myFlags.zmaxval = atof(argv[i + 1]);
    }
    if (!std::strncmp(argv[i], "-vtk", 3)){
      myFlags.vtkOutput = true;
    }

  }
}

void writeUsageHelper(){

  printf("USAGE: multiFrogReader $inputName\n");
  printf("       inputName = stripped name: (e.g. E_FIELD_010.000.binn");
  printf("                                  NOT: E_FIELD_010.000.bin.000)\n");
  printf("       -lockx (to lock the X coord.)\n");
  printf("       -xsample $XSAMPLING (to set the sampling along X)\n");
  printf("       -xmin $XMIN (to set xmin)\n");
  printf("       -xmax $XMAX (to set xmax)\n");
  printf("       -vtk        (to get a vtk output file)\n");
}

void checkFlags(ALL_FLAGS &myFlags, FILE_DATA &fileData, OUTPUT_DATA &outputData){
  myFlags.imaxval[0] = fileData.dataNSize[0];
  myFlags.imaxval[1] = fileData.dataNSize[1];
  myFlags.imaxval[2] = fileData.dataNSize[2];

  if (myFlags.FLAG_xmin)
    myFlags.iminval[0] = findIndexMin(myFlags.xminval, fileData.riCoords[0], fileData.dataNSize[0]);

  if (myFlags.FLAG_xmax)
    myFlags.imaxval[0] = findIndexMax(myFlags.xmaxval, fileData.riCoords[0], fileData.dataNSize[0]);

  if (myFlags.FLAG_ymin)
    myFlags.iminval[1] = findIndexMin(myFlags.yminval, fileData.riCoords[1], fileData.dataNSize[1]);

  if (myFlags.FLAG_ymax)
    myFlags.imaxval[1] = findIndexMax(myFlags.ymaxval, fileData.riCoords[1], fileData.dataNSize[1]);

  if (myFlags.FLAG_zmin)
    myFlags.iminval[2] = findIndexMin(myFlags.zminval, fileData.riCoords[2], fileData.dataNSize[2]);

  if (myFlags.FLAG_zmax)
    myFlags.imaxval[2] = findIndexMax(myFlags.zmaxval, fileData.riCoords[2], fileData.dataNSize[2]);


  for(int c=0; c <3; c++){
    outputData.lockIndex[c]=fileData.dataNSize[c]/2;
    outputData.allocN[c] = (myFlags.imaxval[c] - myFlags.iminval[c])/myFlags.sampling[c];
    if((myFlags.imaxval[c] - myFlags.iminval[c])%myFlags.sampling[c]){
      outputData.allocN[c]++;
    }
  }

  for(int c=0; c<3; c++)
    if(myFlags.FLAG_lockr[c]){
      outputData.lockIndex[c]=findIndexMin(myFlags.lockValue[c], fileData.riCoords[c], fileData.dataNSize[c]);
      //outputData.lockIndex[c]=fileData.dataNSize[c]/2;
      myFlags.iminval[c] = outputData.lockIndex[c];
      myFlags.imaxval[c] = myFlags.iminval[c] + 1;
      outputData.allocN[c]=1;
    }
}

void updateFileIfEOF(std::ifstream &file_bin, int &fileId, std::string strippedName){
  char dummyBuffer[12];
  file_bin.read(dummyBuffer, 12);
  if(file_bin.eof() ){
    file_bin.close();
    fileId ++;
    std::string currentFileName = composeFileName(strippedName,fileId);
    std::cout << "\nApro il file >>> " << currentFileName << std::endl;
    file_bin.open(currentFileName.c_str(), std::ios::binary | std::ios::in);
  }
  else{
    file_bin.seekg(-12, std::ios::cur);
  }
}

void readFromFile(std::ifstream &file_bin, int *buffer, int N, bool doSwap){
  file_bin.read((char*)buffer, N * sizeof(int));
  if(doSwap)
    swap_endian( buffer,N);
}
void readFromFile(std::ifstream &file_bin, float *buffer, int N, bool doSwap){
  file_bin.read((char*)buffer, N * sizeof(float));
  if(doSwap)
    swap_endian( buffer,N);
}

void printVTKFile(OUTPUT_DATA &outputData, ALL_FLAGS &myFlags, FILE_DATA fileData, std::string outputfileName){

  uint64_t totPts = outputData.allocN[0] * outputData.allocN[1]* outputData.allocN[2];
  double dr[3];
  if(!is_big_endian())
    swap_endian(outputData.savedFields, outputData.size);
  dr[0]=fileData.riCoords[0][myFlags.sampling[0]+myFlags.iminval[0]]-fileData.riCoords[0][myFlags.iminval[0]];
  dr[1]=fileData.riCoords[1][myFlags.sampling[1]+myFlags.iminval[1]]-fileData.riCoords[1][myFlags.iminval[1]];
  dr[2]=fileData.riCoords[2][myFlags.sampling[2]+myFlags.iminval[2]]-fileData.riCoords[2][myFlags.iminval[2]];
  for(int c =0; c < 3; c++){
    if(outputData.allocN[c]==1)
      dr[c]=0;
  }
  std::string compNames[fileData.Ncomp];
  if(fileData.Ncomp==1){
    compNames[0] = "scalar";
  }
  else if(fileData.Ncomp==3){
    compNames[0] = "Vx";
    compNames[1] = "Vy";
    compNames[2] = "Vz";
  }


  std::ofstream file_vtk;
  file_vtk.open(outputfileName.c_str());
  std::stringstream bufstream;

  std::cout << "Writing the fields file\n";
  bufstream << "# vtk DataFile Version 2.0\n";
  bufstream << "titolo mio\n";
  bufstream << "BINARY\n";
  bufstream << "DATASET STRUCTURED_POINTS\n";
  bufstream << "DIMENSIONS " << outputData.allocN[0] << "  " << outputData.allocN[1] << "  "  << outputData.allocN[2] << std::endl;
  bufstream << "ORIGIN " << fileData.riCoords[0][myFlags.iminval[0]] << "  "  << fileData.riCoords[1][myFlags.iminval[1]] << "  "  <<  fileData.riCoords[2][myFlags.iminval[2]] << std::endl;
  bufstream << "SPACING " << dr[0] << "  "  << dr[1] << "  "  << dr[2] << std::endl;
  bufstream << "POINT_DATA " << totPts << std::endl;
  std::string bufstring = bufstream.str();
  file_vtk.write(bufstring.c_str(), bufstring.length());

  for (int cc = 0; cc < fileData.Ncomp; cc++){
    std::stringstream ssbuf;
    ssbuf << "SCALARS " << compNames[cc]  <<" float 1\n";
    ssbuf << "LOOKUP_TABLE default\n";
    bufstring = ssbuf.str();
    file_vtk.write(bufstring.c_str(), bufstring.length());
    file_vtk.write((char*)(&outputData.savedFields[cc*totPts]), sizeof(float)*totPts);
  }
  file_vtk.close();
}
