
/*******************************************************************************
This file is part of tools_piccante.

tools_piccante is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
tools_piccante is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with tools_piccante.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/

#include <mpi.h>
#include "Wheader.h"

bool flag_swap=false;
bool flag_mass=false;
double mass;



bool flag_with_filters=false;
struct filter{
  int on_what;
  bool is_there_minval;
  bool is_there_maxval;
  double minval;
  double maxval;
};
std::vector<filter> filterList;

bool flag_inputfile=false;
std::string inputfileName;

bool flag_outputfile=false;
std::string outputfileName;

bool flag_vtk=false;

std::string quantitiesNames[NUM_QUANTITIES] = {"X","Y","Z","Px","Py","Pz","Ptot","Ktot", "phi(X-Y)", "theta(r-Z)"};
bool filter_flags[NUM_QUANTITIES] = {false,false,false,false,false,false,false,false, false};
double min_filter[NUM_QUANTITIES] = {0,0,0,0,0,0,0,0,0};
double max_filter[NUM_QUANTITIES] = {0,0,0,0,0,0,0,0,0};

struct parallelData{
  int myRank;
  int nProc;
  MPI_Comm comm;
};
void initializeMPIWorld(parallelData &world, int *narg, char ***args);

void parseArgs(int narg, char **args, parallelData pdata);
void checkFlagsConsistence(parallelData pdata);
int getAndCheckFileNumber(parallelData pdata, std::string strippedFileName);
int howManyFilesExist(std::string strippedFileName);
long long int getAndCheckFileNumberOfParticles(parallelData pdata, std::string fileName);
long long int howLongIsInputFile(std::string fileName);
void splitCommunicatorFillNewParallelData(parallelData  parent, parallelData &child, int color);
void message(parallelData pdata, std::string msg);
void errorMessage(parallelData pdata, std::string msg);
std::string composeFileName(std::string strippedFileName, int fileId);
long long int calcParticlesToRead(parallelData pdata, long long particleTotalNumber);
void fillComponentsValues(float *components, float *coordinates, float &weight);
void swap_endian_float_array(float* in_f, int n);
void swap_endian_double_array(double* in_f, int n);
void resetMinimaAndMaxima();
void  fillNumberParticlesToRead(parallelData fileGroup, long long *particlesToRead, long long myparticlesToRead);
MPI_Offset findDispForSetView(parallelData fileGroup,long long *particlesToRead);
int findMaxNumberOfReadings(parallelData fileGroup, long long *particlesToRead);
bool shouldIFindExtrems();
std::string printPlotLimits();
void checkMinimaConsistency(parallelData world);
void increasePlotExtremsBy(float factor);

void read_next_chunk(MPI_File myFile, long long numreader, std::string outFileName, parallelData comunicatore);

const int readLength = 10000000;

double mincomponents[NUM_QUANTITIES];
double maxcomponents[NUM_QUANTITIES];

int is_big_endian();
int main(int narg, char **args)
{
  parallelData world;
  initializeMPIWorld(world, &narg, &args);
  parseArgs(narg,args,world);
  if(!flag_with_filters){
    std::cout<<"NO FILTERS!!!"<<std::endl;
    return 1;
  }
  checkFlagsConsistence(world);

  int numberOfFiles = getAndCheckFileNumber(world, inputfileName);
  int fileId = world.myRank%numberOfFiles;

  parallelData fileGroup;
  splitCommunicatorFillNewParallelData(world, fileGroup, fileId);

  std::string dummyString = composeFileName(inputfileName,fileId);
  char * fileName = new char[dummyString.length() + 1];
  std::strcpy(fileName,dummyString.c_str());

  resetMinimaAndMaxima();

  int readings, reminder, maxReadings;
  MPI_Offset disp;
  {
    long long particleTotalNumber = getAndCheckFileNumberOfParticles(fileGroup, fileName);
    long long myparticlesToRead = calcParticlesToRead(fileGroup, particleTotalNumber);

    long long *particlesToRead = new long long[fileGroup.nProc];
    fillNumberParticlesToRead(fileGroup, particlesToRead, myparticlesToRead);

    readings = myparticlesToRead/readLength;
    reminder = myparticlesToRead%readLength;

    disp = findDispForSetView(fileGroup,particlesToRead);

    int maxReadingsFile = findMaxNumberOfReadings(fileGroup,particlesToRead);

    MPI_Allreduce(&maxReadingsFile, &maxReadings, 1, MPI_INT, MPI_MAX, world.comm);

    delete[] particlesToRead;
  }


  MPI_File theFile;
  MPI_File_open(fileGroup.comm, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL , &theFile);
  MPI_File_set_view(theFile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);





  for (int i = 0; i < readings; i++){
    read_next_chunk(theFile,readLength, outputfileName, world);
  }
  for (int i = 0; i < maxReadings-readings; i++){
    read_next_chunk(theFile, 0, outputfileName, world);
  }


  read_next_chunk(theFile, reminder, outputfileName, world);


  MPI_File_close(&theFile);

  MPI_Finalize();
  return 0;
}

void initializeMPIWorld(parallelData &world, int *narg, char ***args){
  MPI_Init(narg, args);
  MPI_Comm_rank(MPI_COMM_WORLD,&world.myRank);
  MPI_Comm_size(MPI_COMM_WORLD,&world.nProc);
  world.comm = MPI_COMM_WORLD;
}

void parseArgs(int nNumberofArgs, char* pszArgs[], parallelData pdata){
  if(nNumberofArgs<2){
    std::stringstream ss;
    ss << "USAGE: " << std::endl;
    ss << "\ttitan -i inputFile -o outputFile -1D (or) -2D (or) -3D" << std::endl;
    ss <<"\t-first $FIRST_COMP -second $SECOND_COMP -third $THIRD_COMP" << std::endl;
    for(int c=0; c < NUM_QUANTITIES; c ++) {
      ss << c << " = "  << quantitiesNames[c].c_str() << std::endl;
    }
    ss << std::endl << "\t-1min $FIRST_COMP_MIN -1max $FIRST_COMP_MAX" << std::endl;
    ss << "\t -1nbin $FIRST_COMP_NBIN" << std::endl;
    ss << "\t-filter $FILTER_COMP:$MIN:$MAX" << std::endl;
    message(pdata, ss.str());
  }
  for (int i = 1; i < nNumberofArgs; i++){
    if (std::string(pszArgs[i]) == "-swap"){
      flag_swap=true;
    }



    if (std::string(pszArgs[i]) == "-filter"){
      flag_with_filters=true;
      if (i + 1 != nNumberofArgs){
        filter new_filter;

        char split_char = ':';
        std::string bstring;
        bstring = std::string(pszArgs[i+1]);
        std::istringstream split(bstring);
        std::vector<std::string> token;
        for (std::string each; std::getline(split, each, split_char); token.push_back(each));

        if(token.size()!=3){
          errorMessage(pdata, "wrong FILTER provided!");
        }


        new_filter.on_what = atoi(token[0].c_str());


        if (token[1]=="*"){
          new_filter.is_there_minval=false;
        }
        else{
          new_filter.is_there_minval=true;
          new_filter.minval=atof(token[1].c_str());
        }



        if (token[2]=="*"){
          new_filter.is_there_maxval=false;
        }
        else{
          new_filter.is_there_maxval=true;
          new_filter.maxval=atof(token[2].c_str());
        }


        filterList.push_back(new_filter);

        i++;
      }
      else{
        errorMessage(pdata, "FILTER not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-i"){
      flag_inputfile=true;
      if (i + 1 != nNumberofArgs){
        inputfileName = std::string(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "INPUTFILE not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-o"){
      flag_outputfile=true;
      if (i + 1 != nNumberofArgs){
        outputfileName = std::string(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "OUTPUTFILE not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-m"){
      flag_mass=true;
      if (i + 1 != nNumberofArgs){
        mass = atof(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "MASS not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-vtk"){
      flag_vtk=true;
    }
  }
}

void checkFlagsConsistence(parallelData pdata){
  std::stringstream ss;
  if(flag_swap){
    ss << "Endianess swap: ON" << std::endl;
  }
  else{
    ss << "Endianess swap: OFF" << std::endl;
  }


  if(flag_inputfile && inputfileName!=""){
    ss << "Input file: " << inputfileName << std::endl;
  }
  else{
    std::cout << "Input file not provided!" << std::endl;
    exit(1);
  }

  if(flag_outputfile && outputfileName!=""){
    ss << "Output file: " << outputfileName << std::endl;
  }
  else{
    std::cout << "Output file not provided!" << std::endl;
    exit(1);
  }

  if(flag_with_filters){
    for (std::vector<filter>::iterator it = filterList.begin() ; it != filterList.end(); ++it){
      filter fil = *it;
      if(filter_flags[fil.on_what]){

        std::cout << "Error! Multiple filter !" << std::endl;
        exit(1);
      }
      filter_flags[fil.on_what]=true;

      max_filter[fil.on_what]=(fil.is_there_maxval)?(fil.maxval):(VERY_BIG_POS_NUM);
      min_filter[fil.on_what]=(fil.is_there_minval)?(fil.minval):(VERY_BIG_NEG_NUM);

      if(max_filter[fil.on_what]<min_filter[fil.on_what]){
        errorMessage(pdata, "Error in filter definition!");
      }

      ss << "Filter on " << quantitiesNames[fil.on_what] << " : " <<  min_filter[fil.on_what] << " -- " << max_filter[fil.on_what];
      ss << std::endl;
    }
  }
  message(pdata,ss.str());
}
long long int getAndCheckFileNumberOfParticles(parallelData pdata, std::string fileName){
  long long fileLengthInBytes=howLongIsInputFile(fileName);
  if(fileLengthInBytes < 0){
    errorMessage(pdata, fileName + " not found!");
  }
  long long particleTotalNumber = fileLengthInBytes/(sizeof(float)*NUM_COMPONENTS);

  std::stringstream ss;
  ss << "There are: " << particleTotalNumber <<" particles in file "<<fileName<< ".";
  message(pdata, ss.str());

  return particleTotalNumber;
}

long long howLongIsInputFile(std::string fileName){
  std::ifstream file( fileName.c_str(), std::ios::binary | std::ios::ate);
  long long length = file.tellg();
  file.close();
  return length;
}

void splitCommunicatorFillNewParallelData(parallelData  parent, parallelData &child, int color){
  MPI_Comm_split(parent.comm, color, 0, &child.comm);
  MPI_Comm_rank(child.comm,&child.myRank);
  MPI_Comm_size(child.comm,&child.nProc);
}

void swap_endian_float_array(float* in_f, int n){
  if(!flag_swap)
    return;
  int i;
  union {int irep; float frep; char arr[4];}x;
  char buff;
  for(i=0;i<n;i++)
  {
    x.frep=in_f[i];
    buff=x.arr[0];
    x.arr[0]=x.arr[3];
    x.arr[3]=buff;
    buff=x.arr[1];
    x.arr[1]=x.arr[2];
    x.arr[2]=buff;
    in_f[i]=x.frep;
  }
}

void swap_endian_double_array(double* in_f, int n){
  if(is_big_endian())
    return;
  int i;
  union {double frep; char arr[8];}x;
  char buff;
  for(i=0;i<n;i++)
  {
    x.frep=in_f[i];
    buff=x.arr[0];
    x.arr[0]=x.arr[7];
    x.arr[7]=buff;

    buff=x.arr[1];
    x.arr[1]=x.arr[6];
    x.arr[6]=buff;

    buff=x.arr[2];
    x.arr[2]=x.arr[5];
    x.arr[5]=buff;

    buff=x.arr[3];
    x.arr[3]=x.arr[4];
    x.arr[4]=buff;

    in_f[i]=x.frep;
  }
}


void read_next_chunk(MPI_File myFile, long long numreader, std::string outFileName, parallelData comunicatore){
  float* buffer = new float[numreader*NUM_COMPONENTS];
  float* outBuffer = new float[numreader*NUM_COMPONENTS];
  bool isOutside = false;
  if (numreader < 0)
    numreader = 0;
  MPI_File_read_all(myFile, buffer, numreader*NUM_COMPONENTS, MPI_FLOAT, MPI_STATUS_IGNORE);

  //  if(numreader == 0)
  //  return;

  swap_endian_float_array(buffer,NUM_COMPONENTS*numreader);

  float components[NUM_QUANTITIES];
  float weight;


  int numberOfParticlesIn=0;
  for(long long i = 0; i < numreader; i++){
    fillComponentsValues(components, &buffer[NUM_COMPONENTS*i+0], weight);
    isOutside =  false;

    for(int icomp=0; icomp<NUM_QUANTITIES; icomp++){
      if((filter_flags[icomp]&&(components[icomp]<min_filter[icomp] || components[icomp]>max_filter[icomp]))){
        isOutside = true;
      }
    }
    if(!isOutside){
      memcpy((void*)&outBuffer[NUM_COMPONENTS*numberOfParticlesIn+0], (void*)&buffer[NUM_COMPONENTS*i+0],NUM_COMPONENTS*sizeof(float));
      numberOfParticlesIn++;
    }

  }
delete[] buffer;


  int *bufferSize=new int[comunicatore.nProc];
  int myBufferSize = numberOfParticlesIn*NUM_COMPONENTS;
  int maxBufferSize;
  bufferSize[comunicatore.myRank] = myBufferSize;
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, bufferSize, 1, MPI_INT, comunicatore.comm);
  MPI_Allreduce(&myBufferSize, &maxBufferSize, 1, MPI_INT, MPI_MAX, comunicatore.comm);

  if(comunicatore.myRank!=0){
    MPI_Send(outBuffer, myBufferSize, MPI_FLOAT, 0, 11, comunicatore.comm);
    delete[] outBuffer;
  }
  else{
    std::ofstream outfile;
    outfile.open(outFileName.c_str(), std::ofstream::app);
    outfile.write((char*)outBuffer,sizeof(float)*myBufferSize);
    delete[] outBuffer;
    float* recvBuffer = new float[maxBufferSize];
    MPI_Status status;
    for(int procID=1; procID<comunicatore.nProc; procID++){
      MPI_Recv(recvBuffer, bufferSize[procID], MPI_FLOAT, procID, 11, comunicatore.comm, &status);
      outfile.write((char*)recvBuffer,sizeof(float)*bufferSize[procID]);
    }
    delete[] recvBuffer;
    outfile.close();
  }


}

int getAndCheckFileNumber(parallelData pdata, std::string strippedFileName){
  int fileNumber = howManyFilesExist(strippedFileName);
  if(fileNumber <= 0)
    errorMessage(pdata, "Input file not found");

  if(fileNumber>pdata.nProc){
    std::stringstream ss;
    ss << "Too many files! " << "( " << fileNumber << " > " << pdata.nProc << " )";
    errorMessage(pdata, ss.str());
  }
  return fileNumber;
}

int howManyFilesExist(std::string strippedFileName){
  int counter = 0;
  int fID = 0;
  bool hMFEflag = true;

  while(hMFEflag){
    std::string bstring = composeFileName(strippedFileName, fID);
    std::ifstream infile;
    infile.open(bstring.c_str());
    if(infile.good()){
      counter++;
      fID++;
      infile.close();
    }
    else{
      hMFEflag=false;
    }
  }

  return counter;
}


void message(parallelData pdata, std::string msg){
  if (pdata.myRank == 0)
    std::cout << msg << std::endl;
  else
    return;
}

void errorMessage(parallelData pdata, std::string msg){
  if (pdata.myRank == 0){
    std::cout << "ERROR: "  << msg << std::endl;
    std::cout.flush();
  }
  MPI_Barrier(pdata.comm);
  MPI_Finalize();
  exit(1);
}

std::string composeFileName(std::string strippedFileName, int fileId){
  const int numZeroes = 5;
  std::stringstream ss;
  ss << strippedFileName <<"."<< std::setfill('0') << std::setw(numZeroes) << fileId;

  return ss.str();
}

long long int calcParticlesToRead(parallelData pdata, long long particleTotalNumber){
  long long numPart = particleTotalNumber/pdata.nProc;
  long long rem = particleTotalNumber%pdata.nProc;

  if(pdata.myRank < rem)
    numPart++;

  return numPart;
}

void fillComponentsValues(float *components, float *coordinates, float &weight){
  components[0]=coordinates[0];//x
  components[1]=coordinates[1];//y
  components[2]=coordinates[2];//z
  components[3]=coordinates[3];//px
  components[4]=coordinates[4];//py
  components[5]=coordinates[5];//pz
  components[6]=sqrt(components[3]*components[3]+components[4]*components[4]+components[5]*components[5]);//ptot
  components[7]=mass*(sqrt(1.0+components[6]*components[6])-1);
  components[8]=atan2(components[4],components[3])/M_PI*180;
  double rr = sqrt(components[3]*components[3]+components[4]*components[4]);
  components[9]=atan2(components[5],rr)/M_PI*180;
  weight = coordinates[6];
}

void resetMinimaAndMaxima(){
  for (int i = 0; i < NUM_QUANTITIES; i++){
    mincomponents[i]=VERY_BIG_POS_NUM;
    maxcomponents[i]=VERY_BIG_NEG_NUM;
  }
}
void  fillNumberParticlesToRead(parallelData fileGroup, long long *particlesToRead, long long myparticlesToRead){
  particlesToRead[fileGroup.myRank] = myparticlesToRead;
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_LONG_LONG_INT, particlesToRead, 1, MPI_LONG_LONG_INT, fileGroup.comm);
}

MPI_Offset findDispForSetView(parallelData fileGroup,long long *particlesToRead){
  MPI_Offset disp = 0;
  for (int i = 0; i < fileGroup.myRank; i++){
    disp+=particlesToRead[i]*NUM_COMPONENTS*sizeof(float);
  }
  return disp;
}

int findMaxNumberOfReadings(parallelData fileGroup, long long *particlesToRead){
  int maxReadings = 0;
  for (int i = 0; i < fileGroup.nProc; i++){
    int ireadings =particlesToRead[i]/readLength;
    if(ireadings>maxReadings)
      maxReadings=ireadings;
  }
  return maxReadings;
}




int is_big_endian(){
  union {
    int i;
    char c[4];
  } bint = { 0x01020304 };

  return bint.c[0] == 1;
}
