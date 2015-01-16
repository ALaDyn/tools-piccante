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


#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<string>
#include<sstream>
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
long first_bins;

bool flag_second_min=false;
double second_min;
bool flag_second_max=false;
double second_max;
bool flag_second_bins=false;
long second_bins;


bool flag_third_min=false;
double third_min;
bool flag_third_max=false;
double third_max;
bool flag_third_bins=false;
long third_bins;

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

#define NUM_QUANTITIES 9
#define VERY_BIG_POS_NUM +1.0e30
#define VERY_BIG_NEG_NUM -1.0e30
std::string quantitiesNames[NUM_QUANTITIES] = {"X","Y","Z","Px","Py","Pz","Ptot","Ktot", "theta2D"};
bool filter_flags[NUM_QUANTITIES] = {false,false,false,false,false,false,false,false};
double min_filter[NUM_QUANTITIES] = {0,0,0,0,0,0,0,0};
double max_filter[NUM_QUANTITIES] = {0,0,0,0,0,0,0,0};

void parseArgs(int narg, char **args);
void checkFlagsConsistence();
long howLongIsInputFile(std::string fileName);
void drawLoadBar(long i, long Ntot, int sizeBar);

void swap_endian_float_array(float* in_f, int n);

void read_next_extremes(std::ifstream& myFile,long numreader);
void read_next_plot(std::ifstream& myFile, long numreader, double* plotData);

const int readLength = 1000000;

double mincomponents[NUM_QUANTITIES];
double maxcomponents[NUM_QUANTITIES];

int main(int narg, char **args)
{
  long fileLengthInBytes, particleTotalNumber;
  parseArgs(narg,args);
  checkFlagsConsistence();

  fileLengthInBytes = howLongIsInputFile(inputfileName);
  if(fileLengthInBytes <= 0){
    std::cout << "ERROR: file empty or not found" << std::endl;
    exit(1);
  }
  particleTotalNumber = fileLengthInBytes/(sizeof(float)*7);
  std::cout << "There are: "<< particleTotalNumber << " particles in the file." << std::endl;

  if(
     !(
       (flag_1D && flag_first_min && flag_first_max)||
       (flag_2D && flag_first_min && flag_first_max && flag_second_min && flag_second_max)||
       (flag_3D && flag_first_min && flag_first_max && flag_second_min && flag_second_max && flag_third_min && flag_third_max)
       )
     ){
    std::cout << "A preliminary reading is needed to find unspecified plot extremes" << std::endl;
    long readings = particleTotalNumber/readLength;
    long reminder = particleTotalNumber-readings*readLength;
    drawLoadBar(0, particleTotalNumber, 30);
    std::ifstream myFile (inputfileName.c_str(), std::ios::in | std::ios::binary);

    for (int i = 0; i < NUM_QUANTITIES; i++){
      mincomponents[i]=VERY_BIG_POS_NUM;
      maxcomponents[i]=VERY_BIG_NEG_NUM;
    }

    for (long i = 0; i < readings; i++){
      read_next_extremes(myFile,readLength);
      drawLoadBar(i*readLength, particleTotalNumber, 30);
    }
    read_next_extremes(myFile,reminder);
    drawLoadBar(particleTotalNumber, particleTotalNumber, 30);
    std::cout << std::endl;
    myFile.close();
    std::cout << "End of the preliminary reading" << std::endl;
    std::cout<<"NEW LIMITS:"<<std::endl;

    if (flag_1D || flag_2D || flag_3D){
      if(!flag_first_min)first_min=mincomponents[what_first];
      if(!flag_first_max)first_max=maxcomponents[what_first];
      std::cout << "First quantity plot limits: " << first_min << " -- " << first_max << std::endl;
    }
    if (flag_2D || flag_3D){
      if(!flag_second_min)second_min=mincomponents[what_second];
      if(!flag_second_max)second_max=maxcomponents[what_second];
      std::cout << "Second quantity plot limits: " << second_min << " -- " << second_max << std::endl;
    }
    if (flag_3D){
      if(!flag_third_min)third_min=mincomponents[what_third];
      if(!flag_third_max)third_max=maxcomponents[what_third];
      std::cout << "Third quantity plot limits: " << third_min << " -- " << third_max << std::endl;
    }
    std::cout << "Plot limits will be increased by a 5%" << std::endl;

    first_min-=0.05*first_min;
    first_max+=0.05*first_max;

    second_min-=0.05*second_min;
    second_max+=0.05*second_max;

    third_min-=0.05*third_min;
    third_max+=0.05*third_max;


  }

  if(first_min>=first_max){
    std::cout << "ERROR! First quantity limits error"<<std::endl;
    exit(1);
  }
  if((flag_2D||flag_3D)&&(second_min>=second_max)){
    std::cout << "ERROR! Second quantity limits error"<<std::endl;
    exit(1);
  }
  if(flag_3D&&(third_min>=third_max)){
    std::cout << "ERROR! Third quantity limits error"<<std::endl;
    exit(2);
  }

  std::cout<<"Preparing plot data..."<<std::endl;

  if(flag_1D){
    second_bins = 1;
    third_bins = 1;
  }
  else if(flag_2D){
    third_bins = 1;
  }

  double* plotData = new double[first_bins*second_bins*third_bins];

  for(long i = 0; i < first_bins; i++)
    for(long j = 0; j < second_bins; j++)
      for(long k = 0; k < third_bins; k++)
      {
        plotData[i+j*first_bins+k*first_bins*second_bins] = 0.0;
      }

  long readings = particleTotalNumber/readLength;
  long reminder = particleTotalNumber-readings*readLength;
  drawLoadBar(0, particleTotalNumber, 30);
  std::ifstream myFile (inputfileName.c_str(), std::ios::in | std::ios::binary);
  for (long i = 0; i < readings; i++){
    read_next_plot(myFile,readLength,plotData);
    drawLoadBar(i*readLength, particleTotalNumber, 30);
  }
  read_next_plot(myFile,reminder,plotData);
  drawLoadBar(particleTotalNumber, particleTotalNumber, 30);
  std::cout << std::endl;
  myFile.close();
  std::cout<<"Writing plot data on disk..."<<std::endl;

  std::ofstream outfile;


  if(flag_1D){
    outfile.open(outputfileName.c_str());

    long j = 0;
    long k = 0;
    for (long i = 0; i < first_bins; i++){
      double fcoord = (i + i + 1)*0.5/first_bins*(first_max-first_min)+first_min;
      outfile << fcoord << " " << plotData[i+j*first_bins+k*first_bins*second_bins] << std::endl;
    }
    outfile.close();
  }
  else if(flag_2D){
    outfile.open(outputfileName.c_str());

    long k = 0;
    for (long j = 0; j < second_bins; j++)
      for (long i = 0; i < first_bins; i++){
        double fcoord = (i + i + 1)*0.5/first_bins*(first_max-first_min)+first_min;
        double scoord = (j + j + 1)*0.5/second_bins*(second_max-second_min)+second_min;
        outfile << fcoord << " " << scoord << " "<< plotData[i+j*first_bins+k*first_bins*second_bins] << std::endl;
      }
    outfile.close();
  }
  else if(flag_3D)
  {
    if(!flag_vtk)
    {
      outfile.open(outputfileName.c_str());

      for (long k = 0; k < third_bins; k++)
        for (long j = 0; j < second_bins; j++)
          for (long i = 0; i < first_bins; i++){
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
      fprintf(clean_fields, "DIMENSIONS %lu %lu %lu\n", first_bins,second_bins,third_bins  );
      fprintf(clean_fields, "ORIGIN %f %f %f\n", first_min, second_min, third_min);
      double dx  = (first_max-first_min)/first_bins;
      double dy = (second_max-second_min)/second_bins;
      double dz  = (third_max-third_min)/third_bins;

      fprintf(clean_fields, "SPACING %f %f %f\n", dx, dy, dz);
      fprintf(clean_fields, "POINT_DATA %lu\n", first_bins*second_bins*third_bins);
      fprintf(clean_fields, "SCALARS ciccio double 1\n" );
      fprintf(clean_fields, "LOOKUP_TABLE default\n");
      fwrite((void*)plotData, sizeof(double), first_bins*second_bins*third_bins, clean_fields);
      fclose(clean_fields);
    }
  }



  delete[] plotData;

  return 0;
}

void parseArgs(int nNumberofArgs, char* pszArgs[]){
  if(nNumberofArgs<2){
    printf("USAGE:\n");
    printf("\ttitan -i inputFile -o outputFile -1D (or) -2D (or) -3D \n");
    printf("\t-first $FIRST_COMP -second $SECOND_COMP -third $THIRD_COMP\n");
    for(int c=0; c < NUM_QUANTITIES; c ++) {
      printf("%i=%s  ", c, quantitiesNames[c].c_str());
    }
    printf("\n\t-1min $FIRST_COMP_MIN -1max $FIRST_COMP_MAX\n");
    printf("\t -1nbin $FIRST_COMP_NBIN\n");
    printf("\t-filter $FILTER_COMP:$MIN:$MAX\n");
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
        std::cout << "ERROR: FIRST quantity not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-second"){
      flag_second=true;
      if (i + 1 != nNumberofArgs){
        what_second = atoi(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: SECOND quantity not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-third"){
      flag_third=true;
      if (i + 1 != nNumberofArgs){
        what_third = atoi(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: THIRD quantity not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-1min"){
      flag_first_min=true;
      if (i + 1 != nNumberofArgs){
        first_min = atof(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: FIRST min not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-2min"){
      flag_second_min=true;
      if (i + 1 != nNumberofArgs){
        second_min = atof(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: SECOND min not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-3min"){
      flag_third_min=true;
      if (i + 1 != nNumberofArgs){
        third_min = atof(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: THIRD min not provided!" << std::endl;
        exit(1);
      }
    }



    if (std::string(pszArgs[i]) == "-1max"){
      flag_first_max=true;
      if (i + 1 != nNumberofArgs){
        first_max = atof(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: FIRST max not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-2max"){
      flag_second_max=true;
      if (i + 1 != nNumberofArgs){
        second_max = atof(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: SECOND max not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-3max"){
      flag_third_max=true;
      if (i + 1 != nNumberofArgs){
        third_max = atof(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: THIRD max not provided!" << std::endl;
        exit(1);
      }
    }


    if (std::string(pszArgs[i]) == "-1nbin"){
      flag_first_bins=true;
      if (i + 1 != nNumberofArgs){
        first_bins = atoi(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: FIRST bins not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-2nbin"){
      flag_second_bins=true;
      if (i + 1 != nNumberofArgs){
        second_bins = atoi(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: SECOND bins not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-3nbin"){
      flag_third_bins=true;
      if (i + 1 != nNumberofArgs){
        third_bins = atoi(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: THIRD bins not provided!" << std::endl;
        exit(1);
      }
    }


    if (std::string(pszArgs[i]) == "-filter"){
      flag_with_filters=true;
      if (i + 1 != nNumberofArgs){
        filter new_filter;

        char split_char = ':';
        std::istringstream split(std::string(pszArgs[i+1]));
        std::vector<std::string> token;
        for (std::string each; std::getline(split, each, split_char); token.push_back(each));

        if(token.size()!=3){
          std::cout << "ERROR: wrong FILTER provided!" << std::endl;
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
        std::cout << "ERROR: FILTER not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-i"){
      flag_inputfile=true;
      if (i + 1 != nNumberofArgs){
        inputfileName = std::string(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: INPUTFILE not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-o"){
      flag_outputfile=true;
      if (i + 1 != nNumberofArgs){
        outputfileName = std::string(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: OUTPUTFILE not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-m"){
      flag_mass=true;
      if (i + 1 != nNumberofArgs){
        mass = atof(pszArgs[i+1]);
        i++;
      }
      else{
        std::cout << "ERROR: MASS not provided!" << std::endl;
        exit(1);
      }
    }

    if (std::string(pszArgs[i]) == "-vtk"){
      flag_vtk=true;
    }
  }
}

void checkFlagsConsistence(){
  if(flag_swap){
    std::cout << "Endianess swap: ON" << std::endl;
  }
  else{
    std::cout << "Endianess swap: OFF" << std::endl;
  }

  if(flag_1D){
    if(flag_2D || flag_3D){
      std::cout << "ERROR, only one plot dimension can be chosen!" << std::endl;
      exit(1);
    }
    std::cout << "Plot type: 1D" << std::endl;
  }
  else if(flag_2D){
    if(flag_3D){
      std::cout << "ERROR, only one plot dimension can be chosen!" << std::endl;
      exit(1);
    }
    std::cout << "Plot type: 2D" << std::endl;
  }
  else if(flag_3D){
    std::cout << "Plot type: 3D" << std::endl;
  }

  if(flag_1D){
    if(flag_first){

      if(what_first < 0 || what_first >= NUM_QUANTITIES){
        std::cout << "Wrong FIRST quantity setting." << std::endl;
        exit(1);
      }

      std::cout << "First quantity: " << quantitiesNames[what_first] << std::endl;
      std::cout << "Second quantity: --" << std::endl;
      std::cout << "Third quantity: --" << std::endl;

      std::cout << "First quantity plot limits: ";
      if(flag_first_min){
        std::cout << first_min;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout << " -- ";
      if(flag_first_max){
        std::cout << first_max;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout<<std::endl;
      if(!flag_first_bins)
        first_bins=200;
      std::cout << "First quantity bins: " << first_bins;
      if(!flag_first_bins)
        std::cout << " (DEFAULT) ";
      std::cout<<std::endl;

      std::cout << "Second quantity plot limits: --" << std::endl;
      std::cout << "Second quantity bins: --" << std::endl;

      std::cout << "Third quantity plot limits: --" << std::endl;
      std::cout << "Third quantity bins: --" << std::endl;

    }
    else{
      std::cout << "ERROR, FIRST quantity should be provided for 1D plot!" << std::endl;
      exit(1);
    }
  }

  if(flag_2D){
    if(flag_first && flag_second){
      if(what_first < 0 || what_first >= NUM_QUANTITIES){
        std::cout << "Wrong FIRST quantity setting." << std::endl;
        exit(1);
      }
      if(what_second < 0 || what_second >= NUM_QUANTITIES){
        std::cout << "Wrong SECOND quantity setting." << std::endl;
        exit(1);
      }

      std::cout << "First quantity: " << quantitiesNames[what_first] << std::endl;
      std::cout << "Second quantity: " << quantitiesNames[what_second] << std::endl;
      std::cout << "Third quantity: --" << std::endl;

      std::cout << "First quantity plot limits: ";
      if(flag_first_min){
        std::cout << first_min;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout << " -- ";
      if(flag_first_max){
        std::cout << first_max;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout<<std::endl;
      if(!flag_first_bins)
        first_bins=200;
      std::cout << "First quantity bins: " << first_bins;
      if(!flag_first_bins)
        std::cout << " (DEFAULT) ";
      std::cout<<std::endl;

      std::cout << "Second quantity plot limits: ";
      if(flag_second_min){
        std::cout << second_min;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout << " -- ";
      if(flag_second_max){
        std::cout << second_max;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout<<std::endl;
      if(!flag_second_bins)
        second_bins=200;
      std::cout << "Second quantity bins: " << second_bins;
      if(!flag_second_bins)
        std::cout << " (DEFAULT) ";
      std::cout<<std::endl;

      std::cout << "Third quantity plot limits: --" << std::endl;
      std::cout << "Third quantity bins: --" << std::endl;
    }
    else{
      std::cout << "ERROR, FIRST and SECOND quantities should be provided for 2D plot!" << std::endl;
      exit(1);
    }
  }

  if(flag_3D){
    if(flag_first && flag_second && flag_third){

      if(what_first < 0 || what_first >= NUM_QUANTITIES){
        std::cout << "Wrong FIRST quantity setting." << std::endl;
        exit(1);
      }
      if(what_second < 0 || what_second >= NUM_QUANTITIES){
        std::cout << "Wrong SECOND quantity setting." << std::endl;
        exit(1);
      }
      if(what_third < 0 || what_third>= NUM_QUANTITIES){
        std::cout << "Wrong THIRD quantity setting." << std::endl;
        exit(1);
      }

      std::cout << "First quantity: " << quantitiesNames[what_first] << std::endl;
      std::cout << "Second quantity: " << quantitiesNames[what_second] << std::endl;
      std::cout << "Third quantity: " << quantitiesNames[what_third] << std::endl;

      std::cout << "First quantity plot limits: ";
      if(flag_first_min){
        std::cout << first_min;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout << " -- ";
      if(flag_first_max){
        std::cout << first_max;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout<<std::endl;
      if(!flag_first_bins)
        first_bins=200;
      std::cout << "First quantity bins: " << first_bins;
      if(!flag_first_bins)
        std::cout << " (DEFAULT) ";
      std::cout<<std::endl;

      std::cout << "Second quantity plot limits: ";
      if(flag_second_min){
        std::cout << second_min;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout << " -- ";
      if(flag_second_max){
        std::cout << second_max;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout<<std::endl;
      if(!flag_second_bins)
        second_bins=200;
      std::cout << "Second quantity bins: " << second_bins;
      if(!flag_second_bins)
        std::cout << " (DEFAULT) ";
      std::cout<<std::endl;

      std::cout << "Third quantity plot limits: ";
      if(flag_third_min){
        std::cout << third_min;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout << " -- ";
      if(flag_third_max){
        std::cout << third_max;
      }
      else{
        std::cout << "NO LIMIT";
      }
      std::cout<<std::endl;
      if(!flag_third_bins)
        third_bins=200;
      std::cout << "Third quantity bins: " << third_bins;
      if(!flag_third_bins)
        std::cout << " (DEFAULT) ";
      std::cout<<std::endl;
    }
    else{
      std::cout << "ERROR, FIRST, SECOND and THIRD quantities should be provided for 3D plot!" << std::endl;
      exit(1);
    }
  }

  if(flag_inputfile && inputfileName!=""){
    std::cout << "Input file: " << inputfileName << std::endl;
  }
  else{
    std::cout << "Input file not provided!" << std::endl;
    exit(1);
  }

  if(flag_outputfile && outputfileName!=""){
    std::cout << "Output file: " << outputfileName << std::endl;
  }
  else{
    std::cout << "Output file not provided!" << std::endl;
    exit(1);
  }

  if(flag_mass){
    std::cout << "Particle mass: " << mass << std::endl;
  }
  else{
    mass = 1.0;
    std::cout << "Particle mass: " << mass << " (DEFAULT)" << std::endl;
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
        std::cout << "Error in filter definition!" << std::endl;
        exit(1);
      }

      std::cout << "Filter on " << quantitiesNames[fil.on_what] << " : " <<  min_filter[fil.on_what] << " -- " << max_filter[fil.on_what];
      std::cout << std::endl;
    }
  }

}

long howLongIsInputFile(std::string fileName){
  std::ifstream file( fileName.c_str(), std::ios::binary | std::ios::ate);
  long length = file.tellg();
  file.close();
  return length;
}

void swap_endian_float_array(float* in_f, int n)
{
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

void read_next_extremes(std::ifstream& myFile,long numreader){
  char* buffer = new char[numreader*sizeof(float)*7];
  myFile.read(buffer, numreader*sizeof(float)*7);

  if(flag_swap){
    swap_endian_float_array((float*)buffer,7*numreader);
  }

  float* fbuf = (float*)buffer;
  float components[NUM_QUANTITIES];

  for(long i = 0; i < numreader; i++){
    components[0]=fbuf[7*i+0];//x
    components[1]=fbuf[7*i+1];//y
    components[2]=fbuf[7*i+2];//z
    components[3]=fbuf[7*i+3];//px
    components[4]=fbuf[7*i+4];//py
    components[5]=fbuf[7*i+5];//pz
    components[6]=sqrt(components[3]*components[3]+components[4]*components[4]+components[5]*components[5]);//ptot
    components[7]=mass*(sqrt(1.0+components[6]*components[6])-1);//ktot
    components[8]=atan2(components[4],components[3])/M_PI*180;

    for(int j = 0; j < NUM_QUANTITIES; j++){
      if(components[j]>maxcomponents[j]) maxcomponents[j]=components[j];
      if(components[j]<mincomponents[j]) mincomponents[j]=components[j];
    }
  }


  delete[] buffer;
}

void read_next_plot(std::ifstream& myFile, long numreader, double* plotData){
  char* buffer = new char[numreader*sizeof(float)*7];
  myFile.read(buffer, numreader*sizeof(float)*7);

  if(flag_swap){
    swap_endian_float_array((float*)buffer,7*numreader);
  }

  float* fbuf = (float*)buffer;
  float components[NUM_QUANTITIES];
  float weight;

  long ibin,jbin,kbin;

  double first_size = first_max-first_min;
  double second_size = second_max-second_min;
  double third_size = third_max-third_min;

  if(flag_1D){
    for(long i = 0; i < numreader; i++){
      components[0]=fbuf[7*i+0];//x
      components[1]=fbuf[7*i+1];//y
      components[2]=fbuf[7*i+2];//z
      components[3]=fbuf[7*i+3];//px
      components[4]=fbuf[7*i+4];//py
      components[5]=fbuf[7*i+5];//pz
      components[6]=sqrt(components[3]*components[3]+components[4]*components[4]+components[5]*components[5]);//ptot
      components[7]=mass*(sqrt(1.0+components[6]*components[6])-1);//ktot
      components[8]=atan2(components[4],components[3])/M_PI*180;
      weight = fbuf[7*i+6];

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
    for(long i = 0; i < numreader; i++){
      components[0]=fbuf[7*i+0];//x
      components[1]=fbuf[7*i+1];//y
      components[2]=fbuf[7*i+2];//z
      components[3]=fbuf[7*i+3];//px
      components[4]=fbuf[7*i+4];//py
      components[5]=fbuf[7*i+5];//pz
      components[6]=sqrt(components[3]*components[3]+components[4]*components[4]+components[5]*components[5]);//ptot
      components[7]=mass*(sqrt(1.0+components[6]*components[6])-1);//ktot
      components[8]=atan2(components[4],components[3])/M_PI*180;
      weight = fbuf[7*i+6];


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
    for(long i = 0; i < numreader; i++){
      components[0]=fbuf[7*i+0];//x
      components[1]=fbuf[7*i+1];//y
      components[2]=fbuf[7*i+2];//z
      components[3]=fbuf[7*i+3];//px
      components[4]=fbuf[7*i+4];//py
      components[5]=fbuf[7*i+5];//pz
      components[6]=sqrt(components[3]*components[3]+components[4]*components[4]+components[5]*components[5]);//ptot
      components[7]=mass*(sqrt(1.0+components[6]*components[6])-1);//ktot
      components[8]=atan2(components[4],components[3])/M_PI*180;
      weight = fbuf[7*i+6];


      if(flag_with_filters){
        for(int icomp=0; icomp<NUM_QUANTITIES; icomp++){
          if(filter_flags[icomp]&&(components[icomp]<min_filter[icomp] || components[icomp]>max_filter[icomp])){
            weight = 0;
          }
        }
      }

      kbin=(long)(third_bins-1)*(components[what_third]-third_min)/third_size;
      jbin=(long)(second_bins-1)*(components[what_second]-second_min)/second_size;
      ibin =(long)(first_bins-1)*(components[what_first]-first_min)/first_size;
      if (ibin >= 0 && ibin < first_bins && jbin>=0 && jbin < second_bins && kbin >= 0 && kbin < third_bins)
        plotData[ibin+jbin*first_bins+kbin*first_bins*second_bins] += weight;
    }
  }
  delete[] buffer;
}

void drawLoadBar(long i, long Ntot, int sizeBar){

  std::cout << "\r";

  float ratio = i /((float)Ntot);

  int numSymbols = (int) sizeBar*ratio;

  std::cout << std::setw(3) << (int)(ratio*100) << "% [";

  int j;
  for (j = 0; j < numSymbols; j++){
    std::cout << "=";
  }
  for (; j < sizeBar; j++){
    std::cout << " ";
  }

  std::cout <<"]";
  std::cout.flush();
}


