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

#include<mpi.h>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<string>
#include<sstream>
#include<istream>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<cstring>
#include<vector>

bool flag_swap=false;

bool flag_1D=false;
bool flag_2D=false;
bool flag_3D=false;

bool flag_first=false;
int what_first;

bool flag_second=false;
int what_second;

bool flag_third=false;
int what_third;

bool flag_first_min=false;
double first_min;
bool flag_first_max=false;
double first_max;
bool flag_first_bins=false;
long long first_bins;

bool flag_second_min=false;
double second_min;
bool flag_second_max=false;
double second_max;
bool flag_second_bins=false;
long long second_bins;


bool flag_third_min=false;
double third_min;
bool flag_third_max=false;
double third_max;
bool flag_third_bins=false;
long long third_bins;

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

bool flag_mass=false;
double mass;

bool flag_vtk=false;

#define NUM_COMPONENTS 7
#define NUM_QUANTITIES 10
#define VERY_BIG_POS_NUM +1.0e30
#define VERY_BIG_NEG_NUM -1.0e30
std::string quantitiesNames[NUM_QUANTITIES] = {"X","Y","Z","Px","Py","Pz","Ptot","Ktot", "phi(X-Y)", "theta(r-Z)"};
bool filter_flags[NUM_QUANTITIES] = {false,false,false,false,false,false,false,false, false};
double min_filter[NUM_QUANTITIES] = {0,0,0,0,0,0,0,0,0};
double max_filter[NUM_QUANTITIES] = {0,0,0,0,0,0,0,0,0};
#define INCREASE_PLOTEXTREMS_FACTOR 0.05

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
void resetMinimaAndMaxima();
void  fillNumberParticlesToRead(parallelData fileGroup, long long *particlesToRead, long long myparticlesToRead);
MPI_Offset findDispForSetView(parallelData fileGroup,long long *particlesToRead);
int findMaxNumberOfReadings(parallelData fileGroup, long long *particlesToRead);
bool shouldIFindExtrems();
std::string printPlotLimits();
void checkMinimaConsistency(parallelData world);
void increasePlotExtremsBy(float factor);

void read_next_extremes(MPI_File myFile,long long numreader);
void read_next_plot(MPI_File myFile, long long numreader, double* plotData);

const int readLength = 1000000;

double mincomponents[NUM_QUANTITIES];
double maxcomponents[NUM_QUANTITIES];


int main(int narg, char **args)
{ 
  parallelData world;
  initializeMPIWorld(world, &narg, &args);
  parseArgs(narg,args,world);
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

    maxReadings = findMaxNumberOfReadings(fileGroup,particlesToRead);

    delete[] particlesToRead;
  }


  if(shouldIFindExtrems()){
    message(world, "A preliminary reading is needed to find unspecified plot extremes");

    MPI_File theFile;
    MPI_File_open(fileGroup.comm, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL , &theFile);
    MPI_File_set_view(theFile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);

    for (int i = 0; i < readings; i++){
      read_next_extremes(theFile,readLength);
    }
    read_next_extremes(theFile,reminder);
    for (int i = 0;i < maxReadings-readings; i++){
      read_next_extremes(theFile, 0);
    }
    MPI_File_close(&theFile);

    MPI_Allreduce(MPI_IN_PLACE, mincomponents, NUM_QUANTITIES, MPI_DOUBLE, MPI_MIN, world.comm);
    MPI_Allreduce(MPI_IN_PLACE, maxcomponents, NUM_QUANTITIES, MPI_DOUBLE, MPI_MAX, world.comm);

    std::string msgExtrems = printPlotLimits();
    message(world, msgExtrems);
    increasePlotExtremsBy(INCREASE_PLOTEXTREMS_FACTOR);
  }

  checkMinimaConsistency(world);
  message(world, "Preparing plot data...");


  double* plotData = new double[first_bins*second_bins*third_bins];

  for( long long i = 0; i < first_bins; i++)
    for(long long j = 0; j < second_bins; j++)
      for(long long k = 0; k < third_bins; k++)
      {
        plotData[i+j*first_bins+k*first_bins*second_bins] = 0.0;
      }
  MPI_File theFile;
  MPI_File_open(fileGroup.comm, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL , &theFile);
  MPI_File_set_view(theFile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);

  for (int i = 0; i < readings; i++){
    read_next_plot(theFile,readLength,plotData);
  }
  read_next_plot(theFile, reminder,plotData);
  for (int i = 0; i < maxReadings-readings; i++){
    read_next_plot(theFile,readLength,plotData);
  }

  MPI_Allreduce(MPI_IN_PLACE, plotData, first_bins*second_bins*third_bins, MPI_DOUBLE, MPI_SUM, world.comm);

  MPI_File_close(&theFile);
  std::stringstream ss;
  ss.str(std::string());
  ss << std::endl;
  ss<<"Writing plot data on disk..."<<std::endl;
  message(world,ss.str());

  if(world.myRank == 0){

    std::ofstream outfile;

    if(flag_1D){
      outfile.open(outputfileName.c_str());

      long long j = 0;
      long long k = 0;
      for (long long i = 0; i < first_bins; i++){
        double fcoord = (i + i + 1)*0.5/first_bins*(first_max-first_min)+first_min;
        outfile << fcoord << " " << plotData[i+j*first_bins+k*first_bins*second_bins] << std::endl;
      }
      outfile.close();
    }
    else if(flag_2D){
      outfile.open(outputfileName.c_str());

      long long k = 0;
      std::stringstream bufstream;
      for (long long j = 0; j < second_bins; j++){
        for (long long i = 0; i < first_bins; i++){
          double fcoord = (i + i + 1)*0.5/first_bins*(first_max-first_min)+first_min;
          double scoord = (j + j + 1)*0.5/second_bins*(second_max-second_min)+second_min;
          bufstream << fcoord << " " << scoord << " "<< plotData[i+j*first_bins+k*first_bins*second_bins] << std::endl;
        }
      }
      std::string bufstring = bufstream.str();
      outfile.write(bufstring.c_str(), bufstring.length());
      outfile.close();
    }
    else if(flag_3D)
    {
      if(!flag_vtk)
      {
        outfile.open(outputfileName.c_str());

        for (long long k = 0; k < third_bins; k++)
          for (long long j = 0; j < second_bins; j++)
            for (long long i = 0; i < first_bins; i++){
              double fcoord = (i + i + 1)*0.5/first_bins*(first_max-first_min)+first_min;
              double scoord = (j + j + 1)*0.5/second_bins*(second_max-second_min)+second_min;
              double tcoord = (k + k + 1)*0.5/third_bins*(third_max-third_min)+third_min;
              outfile << fcoord << " " << scoord << " "<< tcoord << " " << plotData[i+j*first_bins+k*first_bins*second_bins] << std::endl;
            }
        outfile.close();
      }
      else{
        FILE *clean_fields=fopen(outputfileName.c_str(),"wb");

        fprintf(clean_fields, "# vtk DataFile Version 2.0\n");
        fprintf(clean_fields, "titolo mio\n");
        fprintf(clean_fields, "BINARY\n");
        fprintf(clean_fields, "DATASET STRUCTURED_POINTS\n");
        fprintf(clean_fields, "DIMENSIONS %lld %lld %lld\n", first_bins,second_bins,third_bins);
        fprintf(clean_fields, "ORIGIN %f %f %f\n", first_min, second_min, third_min);
        double dx  = (first_max-first_min)/first_bins;
        double dy = (second_max-second_min)/second_bins;
        double dz  = (third_max-third_min)/third_bins;

        fprintf(clean_fields, "SPACING %f %f %f\n", dx, dy, dz);
        fprintf(clean_fields, "POINT_DATA %lld\n", first_bins*second_bins*third_bins);
        fprintf(clean_fields, "SCALARS titanZYX double 1\n" );
        fprintf(clean_fields, "LOOKUP_TABLE default\n");
        fwrite((void*)plotData, sizeof(double), first_bins*second_bins*third_bins, clean_fields);
        fclose(clean_fields);
      }
    }

  }


  delete[] plotData;
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

    if (std::string(pszArgs[i]) == "-1D"){
      flag_1D=true;
    }

    if (std::string(pszArgs[i]) == "-2D"){
      flag_2D=true;
    }

    if (std::string(pszArgs[i]) == "-3D"){
      flag_3D=true;
    }

    if (std::string(pszArgs[i]) == "-first"){
      flag_first=true;
      if (i + 1 != nNumberofArgs){
        what_first = atoi(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "FIRST quantity not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-second"){
      flag_second=true;
      if (i + 1 != nNumberofArgs){
        what_second = atoi(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "SECOND quantity not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-third"){
      flag_third=true;
      if (i + 1 != nNumberofArgs){
        what_third = atoi(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "THIRD quantity not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-1min"){
      flag_first_min=true;
      if (i + 1 != nNumberofArgs){
        first_min = atof(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "FIRST min not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-2min"){
      flag_second_min=true;
      if (i + 1 != nNumberofArgs){
        second_min = atof(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "SECOND min not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-3min"){
      flag_third_min=true;
      if (i + 1 != nNumberofArgs){
        third_min = atof(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "THIRD min not provided!");
      }
    }



    if (std::string(pszArgs[i]) == "-1max"){
      flag_first_max=true;
      if (i + 1 != nNumberofArgs){
        first_max = atof(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "FIRST max not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-2max"){
      flag_second_max=true;
      if (i + 1 != nNumberofArgs){
        second_max = atof(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "SECOND max not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-3max"){
      flag_third_max=true;
      if (i + 1 != nNumberofArgs){
        third_max = atof(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "THIRD max not provided!");
      }
    }


    if (std::string(pszArgs[i]) == "-1nbin"){
      flag_first_bins=true;
      if (i + 1 != nNumberofArgs){
        first_bins = atoi(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "FIRST bins not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-2nbin"){
      flag_second_bins=true;
      if (i + 1 != nNumberofArgs){
        second_bins = atoi(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "SECOND bins not provided!");
      }
    }

    if (std::string(pszArgs[i]) == "-3nbin"){
      flag_third_bins=true;
      if (i + 1 != nNumberofArgs){
        third_bins = atoi(pszArgs[i+1]);
        i++;
      }
      else{
        errorMessage(pdata, "THIRD bins not provided!");
      }
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

  if(flag_1D){
    second_bins = 1;
    third_bins = 1;
    if(flag_2D || flag_3D){
      errorMessage(pdata, "Only one plot dimension can be chosen!");
    }
    ss << "Plot type: 1D" << std::endl;
  }
  else if(flag_2D){
    third_bins = 1;
    if(flag_3D){
      errorMessage(pdata, "Only one plot dimension can be chosen!");
    }
    ss << "Plot type: 2D" << std::endl;
  }
  else if(flag_3D){
    ss << "Plot type: 3D" << std::endl;
  }

  if(flag_1D){
    if(flag_first){

      if(what_first < 0 || what_first >= NUM_QUANTITIES){
        errorMessage(pdata, "Wrong FIRST quantity setting.");
      }

      ss << "First quantity: " << quantitiesNames[what_first] << std::endl;
      ss << "Second quantity: --" << std::endl;
      ss << "Third quantity: --" << std::endl;

      ss << "First quantity plot limits: ";
      if(flag_first_min){
        ss << first_min;
      }
      else{
        ss << "NO LIMIT";
      }
      ss << " -- ";
      if(flag_first_max){
        ss << first_max;
      }
      else{
        ss << "NO LIMIT";
      }
      ss<<std::endl;
      if(!flag_first_bins)
        first_bins=200;
      ss << "First quantity bins: " << first_bins;
      if(!flag_first_bins)
        ss << " (DEFAULT) ";
      ss<<std::endl;

      ss << "Second quantity plot limits: --" << std::endl;
      ss << "Second quantity bins: --" << std::endl;

      ss << "Third quantity plot limits: --" << std::endl;
      ss << "Third quantity bins: --" << std::endl;

    }
    else{
      errorMessage(pdata, "FIRST quantity should be provided for 1D plot!");
    }
  }

  if(flag_2D){
    if(flag_first && flag_second){
      if(what_first < 0 || what_first >= NUM_QUANTITIES){
        errorMessage(pdata, "Wrong FIRST quantity setting.");
      }
      if(what_second < 0 || what_second >= NUM_QUANTITIES){
        errorMessage(pdata, "Wrong SECOND quantity setting.");
      }

      ss << "First quantity: " << quantitiesNames[what_first] << std::endl;
      ss << "Second quantity: " << quantitiesNames[what_second] << std::endl;
      ss << "Third quantity: --" << std::endl;

      ss << "First quantity plot limits: ";
      if(flag_first_min){
        ss << first_min;
      }
      else{
        ss << "NO LIMIT";
      }
      ss << " -- ";
      if(flag_first_max){
        ss << first_max;
      }
      else{
        ss << "NO LIMIT";
      }
      ss<<std::endl;
      if(!flag_first_bins)
        first_bins=200;
      ss << "First quantity bins: " << first_bins;
      if(!flag_first_bins)
        ss << " (DEFAULT) ";
      ss<<std::endl;

      ss << "Second quantity plot limits: ";
      if(flag_second_min){
        ss << second_min;
      }
      else{
        ss << "NO LIMIT";
      }
      ss << " -- ";
      if(flag_second_max){
        ss << second_max;
      }
      else{
        ss << "NO LIMIT";
      }
      ss<<std::endl;
      if(!flag_second_bins)
        second_bins=200;
      ss << "Second quantity bins: " << second_bins;
      if(!flag_second_bins)
        ss << " (DEFAULT) ";
      ss<<std::endl;

      ss << "Third quantity plot limits: --" << std::endl;
      ss << "Third quantity bins: --" << std::endl;
    }
    else{
      errorMessage(pdata, "FIRST and SECOND quantities should be provided for 2D plot!");
    }
  }

  if(flag_3D){
    if(flag_first && flag_second && flag_third){

      if(what_first < 0 || what_first >= NUM_QUANTITIES){
        errorMessage(pdata, "Wrong FIRST quantity setting.");
      }
      if(what_second < 0 || what_second >= NUM_QUANTITIES){
        errorMessage(pdata, "Wrong SECOND quantity setting.");
      }
      if(what_third < 0 || what_third>= NUM_QUANTITIES){
        errorMessage(pdata, "Wrong THIRD quantity setting.");
      }

      ss << "First quantity: " << quantitiesNames[what_first] << std::endl;
      ss << "Second quantity: " << quantitiesNames[what_second] << std::endl;
      ss << "Third quantity: " << quantitiesNames[what_third] << std::endl;

      ss << "First quantity plot limits: ";
      if(flag_first_min){
        ss << first_min;
      }
      else{
        ss << "NO LIMIT";
      }
      ss << " -- ";
      if(flag_first_max){
        ss << first_max;
      }
      else{
        ss << "NO LIMIT";
      }
      ss<<std::endl;
      if(!flag_first_bins)
        first_bins=200;
      ss << "First quantity bins: " << first_bins;
      if(!flag_first_bins)
        ss << " (DEFAULT) ";
      ss<<std::endl;

      ss << "Second quantity plot limits: ";
      if(flag_second_min){
        ss << second_min;
      }
      else{
        ss << "NO LIMIT";
      }
      ss << " -- ";
      if(flag_second_max){
        ss << second_max;
      }
      else{
        ss << "NO LIMIT";
      }
      ss<<std::endl;
      if(!flag_second_bins)
        second_bins=200;
      ss << "Second quantity bins: " << second_bins;
      if(!flag_second_bins)
        ss << " (DEFAULT) ";
      ss<<std::endl;

      ss << "Third quantity plot limits: ";
      if(flag_third_min){
        ss << third_min;
      }
      else{
        ss << "NO LIMIT";
      }
      ss << " -- ";
      if(flag_third_max){
        ss << third_max;
      }
      else{
        ss << "NO LIMIT";
      }
      ss<<std::endl;
      if(!flag_third_bins)
        third_bins=200;
      ss << "Third quantity bins: " << third_bins;
      if(!flag_third_bins)
        ss << " (DEFAULT) ";
      ss<<std::endl;
    }
    else{
      errorMessage(pdata, "FIRST, SECOND and THIRD quantities should be provided for 3D plot!");
    }
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

  if(flag_mass){
    ss << "Particle mass: " << mass << std::endl;
  }
  else{
    mass = 1.0;
    ss << "Particle mass: " << mass << " (DEFAULT)" << std::endl;
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

void swap_endian_float_array(float* in_f, int n)
{
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

void read_next_extremes(MPI_File myFile,long long numreader){
  float* buffer = new float[numreader*NUM_COMPONENTS];

  if (numreader < 0)
    numreader = 0;

  MPI_File_read_all(myFile, (char*)buffer, numreader*NUM_COMPONENTS, MPI_FLOAT, MPI_STATUS_IGNORE);

  if(numreader == 0)
    return;

  swap_endian_float_array(buffer,NUM_COMPONENTS*numreader);

  float components[NUM_QUANTITIES];
  float weight;
  for(long long i = 0; i < numreader; i++){
    fillComponentsValues(components, &buffer[NUM_COMPONENTS*i+0], weight);

    for(int c = 0; c < NUM_QUANTITIES; c++){
      if(components[c]>maxcomponents[c]) maxcomponents[c]=components[c];
      if(components[c]<mincomponents[c]) mincomponents[c]=components[c];
    }
  }
  delete[] buffer;
}

void read_next_plot(MPI_File myFile, long long numreader, double* plotData){
  float* buffer = new float[numreader*NUM_COMPONENTS];

  if (numreader < 0)
    numreader = 0;
  MPI_File_read_all(myFile, buffer, numreader*NUM_COMPONENTS, MPI_FLOAT, MPI_STATUS_IGNORE);

  if(numreader == 0)
    return;

  swap_endian_float_array(buffer,NUM_COMPONENTS*numreader);

  float components[NUM_QUANTITIES];
  float weight;

  long long ibin,jbin,kbin;

  double first_size = first_max-first_min;
  double second_size = second_max-second_min;
  double third_size = third_max-third_min;

  if(flag_1D){
    for(long long i = 0; i < numreader; i++){
      fillComponentsValues(components, &buffer[NUM_COMPONENTS*i+0], weight);

      if(flag_with_filters){
        for(int icomp=0; icomp<NUM_QUANTITIES; icomp++){
          if(filter_flags[icomp]&&(components[icomp]<min_filter[icomp] || components[icomp]>max_filter[icomp])){
            weight = 0;
          }
        }
      }
      kbin=0;
      jbin=0;
      ibin = (first_bins-1)*(components[what_first]-first_min)/first_size;
      if (ibin >= 0 && ibin < first_bins)
        plotData[ibin+jbin*first_bins+kbin*first_bins*second_bins] += weight;
    }
  }
  else if(flag_2D){
    for(long long i = 0; i < numreader; i++){
      fillComponentsValues(components, &buffer[NUM_COMPONENTS*i+0], weight);

      if(flag_with_filters){
        for(int icomp=0; icomp<NUM_QUANTITIES; icomp++){
          if(filter_flags[icomp]&&(components[icomp]<min_filter[icomp] || components[icomp]>max_filter[icomp])){
            weight = 0;
          }
        }
      }
      kbin=0;
      jbin=(second_bins-1)*(components[what_second]-second_min)/second_size;
      ibin = (first_bins-1)*(components[what_first]-first_min)/first_size;
      if (ibin >= 0 && ibin < first_bins && jbin>=0 && jbin < second_bins)
        plotData[ibin+jbin*first_bins+kbin*first_bins*second_bins] += weight;
    }
  }
  else{
    for(long long i = 0; i < numreader; i++){
      fillComponentsValues(components, &buffer[NUM_COMPONENTS*i+0], weight);

      if(flag_with_filters){
        for(int icomp=0; icomp<NUM_QUANTITIES; icomp++){
          if(filter_flags[icomp]&&(components[icomp]<min_filter[icomp] || components[icomp]>max_filter[icomp])){
            weight = 0;
          }
        }
      }

      kbin=(long long)(third_bins-1)*(components[what_third]-third_min)/third_size;
      jbin=(long long)(second_bins-1)*(components[what_second]-second_min)/second_size;
      ibin =(long long)(first_bins-1)*(components[what_first]-first_min)/first_size;
      if (ibin >= 0 && ibin < first_bins && jbin>=0 && jbin < second_bins && kbin >= 0 && kbin < third_bins)
        plotData[ibin+jbin*first_bins+kbin*first_bins*second_bins] += weight;
    }
  }
  delete[] buffer;
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

bool shouldIFindExtrems(){
  return !(
        (flag_1D && flag_first_min && flag_first_max)||
        (flag_2D && flag_first_min && flag_first_max && flag_second_min && flag_second_max)||
        (flag_3D && flag_first_min && flag_first_max && flag_second_min && flag_second_max && flag_third_min && flag_third_max)
        );
}

std::string printPlotLimits(){
  std::stringstream ss;
  ss << std::endl << "End of the preliminary reading" << std::endl;
  ss <<"NEW LIMITS:"<<std::endl;

  if (flag_1D || flag_2D || flag_3D){
    if(!flag_first_min)first_min=mincomponents[what_first];
    if(!flag_first_max)first_max=maxcomponents[what_first];
    ss << "First quantity plot limits: " << first_min << " -- " << first_max << std::endl;
  }
  if (flag_2D || flag_3D){
    if(!flag_second_min)second_min=mincomponents[what_second];
    if(!flag_second_max)second_max=maxcomponents[what_second];
    ss << "Second quantity plot limits: " << second_min << " -- " << second_max << std::endl;
  }
  if (flag_3D){
    if(!flag_third_min)third_min=mincomponents[what_third];
    if(!flag_third_max)third_max=maxcomponents[what_third];
    ss << "Third quantity plot limits: " << third_min << " -- " << third_max << std::endl;
  }
  ss << "Plot limits will be increased by a 5%" << std::endl;
  return ss.str();
}

void checkMinimaConsistency(parallelData world){
  if(first_min>=first_max){
    errorMessage(world, "First quantity limits error");
  }
  if((flag_2D||flag_3D)&&(second_min>=second_max)){
    errorMessage(world, "Second quantity limits error");
  }
  if(flag_3D&&(third_min>=third_max)){
    errorMessage(world, "Third quantity limits error");
  }
}

void increasePlotExtremsBy(float factor){
  first_min-=factor*(first_max - first_min);
  first_max+=factor*(first_max - first_min);

  second_min-=factor*(second_max - second_min);
  second_max+=factor*(second_max - second_min);

  third_min-=factor*(third_max - third_min);
  third_max+=factor*(third_max - third_min);
}
