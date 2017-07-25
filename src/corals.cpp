
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

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <malloc.h>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <vector>
#include "utilities-tools.h"
#include <random>
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/normal_distribution.hpp>
//#include <boost/random/variate_generator.hpp>
#define _DIMENSIONS 2

//typedef boost::uniform_real<> UniformRealDistribution;
//typedef boost::mt19937 MTGenerator;
//typedef boost::variate_generator<MTGenerator&,UniformRealDistribution> Generator;

#define MAX(x,y)	((x)>(y)?(x):(y))
#define PRINT_FREQUENCY 100

#define Z_SIZE 5
#define Y_SIZE 20
#define RADIUS_SIZE 0.05
#define DENSITY_SPACING 0.1
#define TARGET_LENGTH 10

struct PARTICLE{
    int coord[3];
};

struct GRID{
    float radius;
    int Lr[3];
    int xmax;
    double *density;
    int densNpoints;
};
void initializeParticles(PARTICLE *particle, GRID *grid);
void pushParticleSteps(PARTICLE *particle, GRID *grid);
bool isParticleTouching(PARTICLE *particle, std::vector<PARTICLE> *foam, GRID &grid);
bool isTooHigh(PARTICLE *particle, GRID *grid);

int main(int narg, char **args){
  std::vector<PARTICLE> foam;
  std::vector<PARTICLE>::const_iterator iterator;
  GRID grid;
  PARTICLE particle;
  int seed = time(NULL);
  srand(seed);

  grid.xmax = 0;
  grid.Lr[0] = 0;
  grid.Lr[1] = (int) Y_SIZE/RADIUS_SIZE;
  grid.Lr[2] = (int) Z_SIZE/RADIUS_SIZE;

  while(grid.xmax<(TARGET_LENGTH*1.2/RADIUS_SIZE)){

    initializeParticles(&particle, &grid);
    //std::cout << foam.size() << "  " << particle.coord[0] << "  " << particle.coord[1] << "  " << particle.coord[2] << std::endl;
    while(1){
      pushParticleSteps(&particle,&grid);
      if(isParticleTouching(&particle,&foam, grid) ){
        foam.push_back(particle);
        grid.xmax = MAX(grid.xmax,particle.coord[0]);
        break;
      }
      if(isTooHigh(&particle, &grid)){
        initializeParticles(&particle, &grid);
      }

    }

    if(!(foam.size()%PRINT_FREQUENCY)){
      std::cout << " size = " << foam.size() <<  "\n";
    }

  }

  int Npart=foam.size();
  int pointerSize = Npart*(3+1);
  float *spheresCoords = new float[pointerSize];
  for (int p=0; p < Npart; p++){
    for(int i=0; i<3; i++){
      spheresCoords[p*4+i] = foam[p].coord[i]*RADIUS_SIZE;
    }
    spheresCoords[p*4+3] = RADIUS_SIZE;
  }


  std::ofstream of1;
  of1.open("spheres.txt", std::ofstream::out | std::ofstream::trunc);
  of1 << "#  " << foam.size() << std::endl;
  of1 << "#  " << std::endl;
  of1 << "#  0  0  0 "<< std::endl;
  of1 << "#  " << grid.xmax <<  " " << grid.Lr[1] <<  " " << grid.Lr[2] <<  " " << std::endl;
  for (int p=0; p < foam.size(); p++){
    for(int i=0; i<4; i++){
      of1 << spheresCoords[p*4+i] << " ";
    }
    of1 <<  std::endl;
  }
  of1.close();


}
void initializeParticles(PARTICLE *particle, GRID *grid){
  particle->coord[0] = grid->xmax + 10;
  particle->coord[1] = 0;
  particle->coord[2] = 0;
  for(int i=1; i<_DIMENSIONS; i++){
    float buffer = (rand()*1.0/RAND_MAX);
    particle->coord[i] = (int)(grid->Lr[i]*buffer);
  }
}

void pushParticleSteps(PARTICLE *particle, GRID *grid){
  float PXU, PXD, PYU, PYD, PZU, PZD;
  float Ptot;
  if(_DIMENSIONS == 2){
    PXD = 1;
    PXU = 1;
    PYD = PYU = 1;
    Ptot = PXU+PXD+PYU+PYD;
    float buffer = Ptot*(rand()*1.0/RAND_MAX);

    if(buffer<=PXD)
      particle->coord[0]--;
    else if(buffer<=(PXD+PXU))
      particle->coord[0]++;
    else if(buffer<=(PXD+PXU+PYD))
      particle->coord[1] = (particle->coord[1]-1+grid->Lr[1])%grid->Lr[1];
    else
      particle->coord[1] = (particle->coord[1]+1)%grid->Lr[1];

  }
  else if(_DIMENSIONS == 3){
    PXD = 1.1;
    PXU = 1;
    PYD = PYU = 1;
    PZD = PZU = 1;
    Ptot = PXU + PXD + PYU + PYD + PZU + PZD;

    float buffer = Ptot*(rand()*1.0/RAND_MAX);
    if(buffer<=PXD)
      particle->coord[0]--;

    else if(buffer<=(PXD+PXU))
      particle->coord[0]++;

    else if(buffer<=(PXD+PXU+PYD))
      particle->coord[1] = (particle->coord[1]-1+grid->Lr[1])%grid->Lr[1];

    else if(buffer<=(PXD+PXU+PYD+PYU))
      particle->coord[1] = (particle->coord[1]+1)%grid->Lr[1];

    else if(buffer<=(PXD+PXU+PYD+PYU+PZD))
      particle->coord[2] = (particle->coord[2]-1+grid->Lr[2])%grid->Lr[2];

    else
      particle->coord[2] = (particle->coord[2]+1)%grid->Lr[2];

  }
}

bool isParticleTouching(PARTICLE *particle, std::vector<PARTICLE> *foam, GRID &grid){
  int dist;
  if(particle->coord[0] == 0)
    return true;
  for (std::vector<PARTICLE>::reverse_iterator i = foam->rbegin(); i != foam->rend(); ++i){
    dist = 0;
    int dummy;

    dist += abs( ((*i).coord[0] - particle->coord[0]) ) ;
    for(int dim =1; dim < _DIMENSIONS; dim ++){
      dummy = abs( ((*i).coord[dim] - particle->coord[dim]) ) ;
      if(dummy == grid.Lr[dim] -1)
        dummy=1;

      dist += dummy;

    }
    if(dist < 2 )
      return true;

  }
  return false;
}

bool isTooHigh(PARTICLE *particle, GRID *grid){
  return (particle->coord[0] > (grid->xmax +20));
}
