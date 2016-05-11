#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string.h>
#include <cstdint>
#include <cmath>

#define SWAPENDIANESS false

#define MAX(A,B) (A>B)?(A):(B)

const double XSIZE = 20;
const double YSIZE = 20;
const double RADIUS = 0.025;

const double FILLING = 0.01;

int is_big_endian(){
  union {
    uint32_t i;
    char c[4];
  } bint = { 0x01020304 };

  return bint.c[0] == 1;
}



void swap_endian(float* in_f, uint64_t n)
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

void swap_endian(double* in_d, int n){
  if(is_big_endian())
    return;
  int i;
  union {double frep; char arr[8];}x;
  char buff;
  for(i=0;i<n;i++){
    x.frep=in_d[i];
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

    in_d[i]=x.frep;
  }
}


int main(int nargs, char** args){
	srand(time(NULL));

	int numpart = (int) (FILLING*(XSIZE*YSIZE)/M_PI/(RADIUS*RADIUS));

	std::cout << "numpart :" << numpart << std::endl; 
	
	double* plist = new double[2*numpart];
	memset (plist,0,2.0*numpart); 

	for (int i = 0; i < numpart; i++){
		double xpos = rand()*XSIZE/RAND_MAX;
		double ypos = rand()*YSIZE/RAND_MAX;

		for (int j = i-1; j >= 0; j--){
			double xj = plist[2*j+0];
			double yj = plist[2*j+1];
			double r2 = (ypos -yj)*(ypos -yj)+(ypos -yj)*(ypos -yj);
			if(r2 < RADIUS*RADIUS){
				i--;
				break;
			}			
		}
		plist[2*i+0] = xpos;
		plist[2*i+1] = ypos;
		
	}
	
	std::ofstream of;
	of.open("spheres.txt");
	for (int i = 0; i < numpart; i++)
		of << plist[2*i+0] << " " << plist[2*i+1] << " " << 2.0*RADIUS << std::endl;
	
	of.close();

	std::cout << numpart << std::endl;
	
		
	float *spheresCoords = new float[numpart*4];
	for (int i = 0; i < numpart; i++){
		spheresCoords[i*4+0] = (float)plist[2*i+0];
		spheresCoords[i*4+1] = (float)plist[2*i+1];	
		spheresCoords[i*4+2] = 0.0f;	
		spheresCoords[i*4+3] = RADIUS;		
	}

	std::ofstream of1;
	float rMin[3];
	rMin[0] = 0.0;
	rMin[1] = 0.0;
	rMin[2] = 0.0;
	float rMax[3];
	rMax[0] = XSIZE;
	rMax[1] = YSIZE;
	rMax[2] = 0.0;

	int dummycount = numpart;

	float fillingFactor=(float)(numpart*M_PI*RADIUS*RADIUS/XSIZE/YSIZE);	
	std::cout << "effective filling: " << fillingFactor << std::endl;

	if(SWAPENDIANESS){
		swap_endian(&dummycount,1);
		swap_endian(&fillingFactor,1);
		swap_endian(rMin,3);
		swap_endian(rMax,3);
		swap_endian(spheresCoords,numpart*4);
		of1.open("spheres_swap.bin", std::ofstream::out | std::ofstream::trunc);
	}
	else{
		of1.open("spheres.bin", std::ofstream::out | std::ofstream::trunc);	
	}
	
	of1.write((char*)&dummycount,sizeof(int));
	of1.write((char*)&fillingFactor,sizeof(float));
	of1.write((char*)rMin,sizeof(float)*3);
	of1.write((char*)rMax,sizeof(float)*3);
	of1.write((char*)spheresCoords,sizeof(float)*4*numpart);
	of1.close();

	delete[] spheresCoords;
	delete[] plist;
	
	return 0;
}

