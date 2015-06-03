#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string.h>
#include <cstdint>
#include <cmath>

#define SWAPENDIANESS true

#define MAX(A,B) (A>B)?(A):(B)

#define XCELLS 250
#define XCUT 160
#define XTGT 200
#define YCELLS 400
#define ZCELLS 400
#define DIAMETER 0.05

struct particle{
	int i;
	int j;
	int k;
};

bool checkBoundary(particle p,char* mat){
	if (p.i > 0 && mat[p.i-1+p.j*XCELLS+p.k*XCELLS*YCELLS]==1)
		return true;

	if (p.i < XCELLS-1 && mat[p.i+1+p.j*XCELLS+p.k*XCELLS*YCELLS]==1)
		return true;

	if (p.j > 0 && mat[p.i+(p.j-1)*XCELLS+p.k*XCELLS*YCELLS]==1)
		return true;

	if (p.j == 0 && mat[p.i+(YCELLS-1)*XCELLS+p.k*XCELLS*YCELLS]==1)
		return true;

	if (p.j < YCELLS-1 && mat[p.i+(p.j+1)*XCELLS+p.k*XCELLS*YCELLS]==1)
		return true;

	if (p.j == YCELLS-1 && mat[p.i+(0)*XCELLS+p.k*XCELLS*YCELLS]==1)
		return true;

	if (p.k > 0 && mat[p.i+p.j*XCELLS+(p.k-1)*XCELLS*YCELLS]==1)
		return true;

	if (p.k == 0 && mat[p.i+p.j*XCELLS+(ZCELLS-1)*XCELLS*YCELLS]==1)
		return true;

	if (p.k < ZCELLS-1 && mat[p.i+p.j*XCELLS+(p.k+1)*XCELLS*YCELLS]==1)
		return true;

	if (p.k == ZCELLS-1 && mat[p.i+p.j*XCELLS+(0)*XCELLS*YCELLS]==1)
		return true;

	return false;
}

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

	char* mat = new char[XCELLS*YCELLS*ZCELLS];
	memset (mat,0,XCELLS*YCELLS*ZCELLS); 

	int imax = 0;

	while(imax < XTGT)
        {
		particle p;
		p.i = imax + 10;
		p.j = rand()%YCELLS;
		p.k = rand()%ZCELLS;		

		//std::cout << p.j << " " << p.k << std::endl;

		while(1){
	
			double pxd = 1.05;
			double pxu = 1.0;
			double pyd = 1.0;
			double pyu = 1.0;
			double pzd = 1.0;
			double pzu = 1.0;
			double ptot = pxd+pxu+pyd+pyu+pzd+pzu;
			double dice = rand()*ptot/RAND_MAX;

			if (dice <= pxd)
				p.i--;
			else if (dice <= pxd+pxu)
				p.i++;	
			else if (dice <= pxd+pxu+pyd)
				p.j--;	
			else if (dice <= pxd+pxu+pyd+pyu)
				p.j++;	
			else if (dice <= pxd+pxu+pyd+pyu+pzd)
				p.k--;	
			else if (dice <= pxd+pxu+pyd+pyu+pzd+pzu)
				p.k++;	


				
			/*			double dice = rand()%7;

			if (dice <= 1)
				p.i--;
			else if (dice == 2)
				p.i++;	
			else if (dice == 3)
				p.j--;	
			else if (dice == 4)
				p.j++; 
			else if (dice == 5)
				p.k--;	
			else if (dice == 6)
			p.k++;*/

			if (p.i == 0){				
				break;			
			}
			if (p.i >= XCELLS-1){
				p.i = imax + 10;
				p.j = rand()%YCELLS;
				p.k = rand()%ZCELLS;
			}

			if(p.j < 0)
				p.j = YCELLS-1;
			else if(p.j > YCELLS-1)
				p.j = 0;

			if(p.k < 0)
				p.k = ZCELLS-1;
			else if(p.k > ZCELLS-1)
				p.k = 0;
			
	
			if (checkBoundary(p,mat)){
				
				break;			
			}

			

		}


		int oldmax = imax;
		imax = MAX(imax, p.i);
		if (imax > oldmax)
			std::cout << imax << std::endl;
		mat[p.i+p.j*XCELLS+p.k*XCELLS*YCELLS] = 1;
		
		
	
	}

	int counter = 0;

	std::ofstream of;
	of.open("coral.txt");
	for (int i = 0; i < XCELLS; i++){
		for (int j = 0; j < YCELLS; j++){
			for (int k = 0; k < ZCELLS; k++){
				if( mat[i+j*XCELLS+k*XCELLS*YCELLS] == 1){
					of << i*DIAMETER << " " << j*DIAMETER << " " << k*DIAMETER << std::endl;		
					counter++;	
				}	
			}		
		}	
	}	

	std::cout << counter << std::endl;
	of.close();

	int rcounter = 0;	
	float *spheresCoords = new float[counter*4];
	float boxVolume = XCUT*YCELLS*ZCELLS;
	float foamVolume = 0.0;

	for (int i = 0; i < XCUT; i++){
		for (int j = 0; j < YCELLS; j++){
			for (int k = 0; k < ZCELLS; k++){
				if( mat[i+j*XCELLS+k*XCELLS*YCELLS] == 1){
					spheresCoords[rcounter*4 + 0] = i*DIAMETER;
					spheresCoords[rcounter*4 + 1] = j*DIAMETER;
					spheresCoords[rcounter*4 + 2] = k*DIAMETER;
					spheresCoords[rcounter*4 + 3] = DIAMETER*0.5;
					foamVolume += 4.0/3.0*M_PI/8.0;
					rcounter++;
				}
			}		
		}
	}
	float fillingFactor = foamVolume/boxVolume;
	std::cout << fillingFactor << " "<< counter << " " << rcounter << std::endl;
	std::cout.flush();
		
	std::ofstream of1;
	float rMin[3];
	rMin[0] = 0.0;
	rMin[1] = 0.0;
	rMin[2] = 0.0;
	float rMax[3];
	rMax[0] = XCUT*DIAMETER;
	rMax[1] = YCELLS*DIAMETER;
	rMax[2] = ZCELLS*DIAMETER;

	int dummycount = rcounter;
	
	if(SWAPENDIANESS){
		swap_endian(&dummycount,1);
		swap_endian(&fillingFactor,1);
		swap_endian(rMin,3);
		swap_endian(rMax,3);
		swap_endian(spheresCoords,rcounter*4);
		of1.open("spheres_swap.bin", std::ofstream::out | std::ofstream::trunc);
	}
	else{
		of1.open("spheres.bin", std::ofstream::out | std::ofstream::trunc);	
	}
	
	std::cout << fillingFactor << std::endl;
	std::cout.flush();	
		
	of1.write((char*)&dummycount,sizeof(int));
	of1.write((char*)&fillingFactor,sizeof(float));
	of1.write((char*)rMin,sizeof(float)*3);
	of1.write((char*)rMax,sizeof(float)*3);
	of1.write((char*)spheresCoords,sizeof(float)*4*rcounter);
	of1.close();

	
	return 0;
}

