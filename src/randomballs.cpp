
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

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 500
#endif

#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string.h>
#include <cstdint>
#include <cmath>
#include "utilities-tools.h"

#define SWAPENDIANESS false

#define MAX(A,B) (A>B)?(A):(B)

const double XSIZE = 20;
const double YSIZE = 20;
const double RADIUS = 0.025;

const double FILLING = 0.01;



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

