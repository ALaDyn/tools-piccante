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

#define Z_MIN 0
#define Z_MAX 5
#define Y_MIN 0
#define Y_MAX 20
#define STEP_SIZE 0.1
#define RADIUS_SIZE 0.05
#define THRESHOLD 0.001
#define DENSITY_SPACING 0.1
#define TARGET_LENGTH 10

struct PARTICLE{
    float coord[3];
    float radius;
};

struct GRID{
    float radius;
    float rMin[3], rMax[3];
    float step;
    float threshold;
    float spacing;
    double *density;
    int densNpoints;
};

void initializeParticles(PARTICLE *particle, GRID *grid){
  particle->radius = grid->radius;
  particle->coord[0] = grid->rMax[0] + grid->step*10;
  particle->coord[1] = 0;
  particle->coord[2] = 0;
  for(int i=1; i<_DIMENSIONS; i++){
    float buffer = (rand()*1.0/RAND_MAX);
    int Nx = (int)(grid->rMax[i] - grid->rMin[i])/grid->step;
    int xint = (int)(Nx*buffer);
    particle->coord[i] = grid->rMin[i] + xint*grid->step;
  }
}

void pushDownParticle(PARTICLE *particle, GRID *grid){
  particle->coord[0] -= grid->step;
}

void pushSideParticle(PARTICLE *particle, GRID *grid){
  float buffer = (rand()*1.0/RAND_MAX);

  float angle = 2*M_PI*buffer-M_PI;
  if(_DIMENSIONS == 2){
    particle->coord[1] += grid->step*angle/M_PI;
  }
  else if(_DIMENSIONS == 3){
    particle->coord[1] += grid->step*cos(angle);
    particle->coord[2] += grid->step*sin(angle);
  }
}

void pushSideParticleSteps(PARTICLE *particle, GRID *grid){
  float buffer = (rand()*1.0/RAND_MAX);

  float angle = 2*M_PI*buffer-M_PI;
  if(_DIMENSIONS == 2){
    if(buffer<=0.5)
      particle->coord[1] -= grid->step;
    else
      particle->coord[1] += grid->step;

  }
  else if(_DIMENSIONS == 3){
    if(buffer<=0.25)
      particle->coord[1] -= grid->step;
    else if(buffer<=0.5)
      particle->coord[1] += grid->step;
    else if(buffer<=0.75)
      particle->coord[2] -= grid->step;
    else
      particle->coord[2] += grid->step;

  }
}
void pushParticleSteps(PARTICLE *particle, GRID *grid){
  float PXU, PXD, PYU, PYD, PZU, PZD;
  float Ptot;
  if(_DIMENSIONS == 2){
    PXD = 1.1;
    PXU = 1;
    PYD = PYU = 1.;
    Ptot = PXU+PXD+PYU+PYD;
    float buffer = Ptot*(rand()*1.0/RAND_MAX);

    if(buffer<=PXD)
      particle->coord[0] -= grid->step;
    else if(buffer<=(PXD+PXU))
      particle->coord[0] += grid->step;
    else if(buffer<=(PXD+PXU+PYD))
      particle->coord[1] -= grid->step;
    else
      particle->coord[1] += grid->step;

  }
  else if(_DIMENSIONS == 3){
    PXD = 1.1;
    PXU = 1;
    PYD = PYU = 1;
    PZD = PZU = 1;
    Ptot = PXU + PXD + PYU + PYD + PZU + PZD;

    float buffer = Ptot*(rand()*1.0/RAND_MAX);
    if(buffer<=PXD)
      particle->coord[0] -= grid->step;

    else if(buffer<=(PXD+PXU))
      particle->coord[0] += grid->step;

    else if(buffer<=(PXD+PXU+PYD))
      particle->coord[1] -= grid->step;

    else if(buffer<=(PXD+PXU+PYD+PYU))
      particle->coord[1] += grid->step;

    else if(buffer<=(PXD+PXU+PYD+PYU+PZD))
      particle->coord[2] -= grid->step;

    else
      particle->coord[2] += grid->step;

  }
}

void pushParticleVariSteps(PARTICLE *particle, GRID *grid){
  float PXU, PXD, PYU, PYD, PZU, PZD;
  float Ptot;
  if(_DIMENSIONS == 2){
    PXD = 1.1;
    PXU = 1;
    PYD = PYU = 1.2;
    Ptot = PXU+PXD+PYU+PYD;
    float buffer = Ptot*(rand()*1.0/RAND_MAX);
    float length = 1.0*(rand()*1.0/RAND_MAX);

    if(buffer<=PXD)
      particle->coord[0] -= length*grid->step;
    else if(buffer<=(PXD+PXU))
      particle->coord[0] += length*grid->step;
    else if(buffer<=(PXD+PXU+PYD))
      particle->coord[1] -= length*grid->step;
    else
      particle->coord[1] += length*grid->step;

  }
  else if(_DIMENSIONS == 3){
    PXD = 2.1;
    PXU = 1;
    PYD = PYU = 1;
    PZD = PZU = 1;
    Ptot = PXU + PXD + PYU + PYD + PZU + PZD;

    float buffer = Ptot*(rand()*1.0/RAND_MAX);
    float length = 1.0*(rand()*1.0/RAND_MAX);
    if(buffer<=PXD)
      particle->coord[0] -= length*grid->step;
    else if(buffer<=(PXD+PXU))
      particle->coord[0] += length*grid->step;
    else if(buffer<=(PXD+PXU+PYD))
      particle->coord[1] -= length*grid->step;
    else if(buffer<=(PXD+PXU+PYD-PYU))
      particle->coord[1] += length*grid->step;
    else if(buffer<=(PXD+PXU+PYD-PYU+PZD))
      particle->coord[2] -= length*grid->step;
    else
      particle->coord[2] += length*grid->step;

  }
}

void checkBoundaries(PARTICLE *particle, GRID *grid){
  for(int d=1; d<_DIMENSIONS; d++){
    if(particle->coord[d] < grid->rMin[d]){
      particle->coord[d] += grid->rMax[d] - grid->rMin[d];
    }
    if(particle->coord[d] >= grid->rMax[d]){
      particle->coord[d] -= grid->rMax[d] - grid->rMin[d];
    }
  }
}

float distance2(PARTICLE *part1, PARTICLE *part2, GRID &grid){
  float r2=0;
  float dr[3];
  float size[3];
  size[0] = grid.rMax[0] - grid.rMin[0];
  size[1] = grid.rMax[1] - grid.rMin[1];
  size[2] = grid.rMax[2] - grid.rMin[2];

  dr[0] = fabs(part1->coord[0]-part2->coord[0]);
  for(int i=1; i<_DIMENSIONS; i++){
    dr[i] = fabs(part1->coord[i]-part2->coord[i]);
    dr[i] -= ((int)(dr[i]/size[i]) + 0.5)*size[i];
  }
  for(int i=0; i<_DIMENSIONS; i++){
    r2 +=  (dr[i])*(dr[i]);
  }
  return r2;
}

bool isParticleTouching(PARTICLE *particle, std::vector<PARTICLE> *foam, GRID &grid){
  float dist;
  //std::vector<PARTICLE>::const_iterator iterator;
  if(particle->coord[0]<=(particle->radius))
    return true;
  for (std::vector<PARTICLE>::reverse_iterator i = foam->rbegin(); i != foam->rend(); ++i){
    dist = distance2( particle, &(*i),  grid);
    if(dist < ( grid.threshold + particle->radius + (*i).radius) * (grid.threshold + particle->radius + (*i).radius) )
      return true;

  }
  return false;
}

void evaluateDensity(std::vector<PARTICLE> &foam, GRID &grid){

  memset( (void*)grid.density, 0 ,grid.densNpoints*sizeof(double) );
  if(_DIMENSIONS == 2){
    for (int p=0; p < foam.size(); p++){
      int ii = (int)(foam[p].coord[0]/grid.spacing);
      if(ii<0 || ii>=grid.densNpoints){
        exit(11);
      }
      grid.density[ii] += (double)(M_PI*foam[p].radius*foam[p].radius);
   }
    for(int i=0; i<grid.densNpoints; i++){
      grid.density[i]/=(grid.rMax[1]-grid.rMin[1])*grid.spacing;
    }
  }
  else if(_DIMENSIONS == 3){
    for (int p=0; p < foam.size(); p++){
      int ii=foam[p].coord[0]/grid.spacing;
      if(ii<0 || ii>=grid.densNpoints){
        exit(11);
      }
      grid.density[ii] += (double)(4.0/3.0*M_PI*foam[p].radius*foam[p].radius*foam[p].radius);

    }
    for(int i=0; i<grid.densNpoints; i++){
      grid.density[i]/=(grid.rMax[1]-grid.rMin[1])*(grid.rMax[2]-grid.rMin[2])*grid.spacing;
    }
  }
}
void printDensity(int ID, GRID &grid){

  std::stringstream ss;
  ss << "dens_" << ID << ".txt";
  std::ofstream of1;
  of1.open(ss.str().c_str(), std::ofstream::out | std::ofstream::trunc);
  for(int i=0; i<grid.densNpoints; i++){
    float x = grid.spacing*i;
    of1 << x << "  " << grid.density[i] << std::endl;
  }
  of1.close();
}
int main(int narg, char **args){
  std::vector<PARTICLE> foam;
  std::vector<PARTICLE>::const_iterator iterator;
  GRID grid;
  PARTICLE particle;
  int seed = time(NULL);
  srand(seed);
  //UniformRealDistribution distribution(0, 1.0);
  //MTGenerator generator;
  //Generator numberGenerator(generator, distribution);
  //generator.seed(0);
  //grid.myGenerator = numberGenerator;

  grid.rMin[0] = 0.0;
  grid.rMax[0] = 0.0;
  grid.rMin[1] = Y_MIN;
  grid.rMax[1] = Y_MAX;
  grid.rMin[2] = Z_MIN;
  grid.rMax[2] = Z_MAX;

  grid.step = STEP_SIZE;
  grid.radius = RADIUS_SIZE;
  grid.threshold = THRESHOLD;
  grid.spacing = DENSITY_SPACING;

  std::cout << "  grid.step = " << grid.step << "  grid.radius = " << grid.radius << "\n";
  std::cout << "  foam.size()  = " << foam.size() << "    particle.coord[0] = " << particle.coord[0] << std::endl;
  for(int i=0; i<_DIMENSIONS; i++){
    std::cout << "  grid.rMin[" << i << "]  = " << grid.rMin[i] << "  grid.rMax[" << i << "]  = " << grid.rMax[i] << std::endl;
  }
  initializeParticles(&particle, &grid);

  std::cout << "  foam.size()  = " << foam.size() << "    particle.coord[0] = " << particle.coord[0] << std::endl;
  std::cout << "  foam.size()  = " << foam.size() << "    particle.coord[0] = " << particle.coord[0] << std::endl;

  while(grid.rMax[0]<(TARGET_LENGTH*1.2)){
    while(1){
#ifdef OLD
      pushDownParticle(&particle,&grid);
      if(isParticleTouching(&particle,&foam, grid) ){
        foam.push_back(particle);
        grid.rMax[0] = MAX(grid.rMax[0],particle.coord[0]);
        break;
      }

      pushSideParticleSteps(&particle,&grid);
      checkBoundaries(&particle, &grid);
      if(isParticleTouching(&particle,&foam, grid) ){
        foam.push_back(particle);
        grid.rMax[0] = MAX(grid.rMax[0],particle.coord[0]);
        break;
      }
#else
      pushParticleSteps(&particle,&grid);
      //pushParticleVariSteps(&particle,&grid);
      checkBoundaries(&particle, &grid);
      if(isParticleTouching(&particle,&foam, grid) ){
        foam.push_back(particle);
        grid.rMax[0] = MAX(grid.rMax[0],particle.coord[0]);
         initializeParticles(&particle, &grid);
        break;
      }
#endif

    }

    if(!(foam.size()%PRINT_FREQUENCY)){
      //grid.densNpoints = grid.rMax[0]/grid.spacing + 5;
      //grid.density = new double[grid.densNpoints];
      //evaluateDensity(foam, grid);
      //printDensity(foam.size(), grid);
      //delete [] grid.density;
      std::cout << " size = " << foam.size() <<  "\n";
    }
    initializeParticles(&particle, &grid);
  }
std::cout << " "<< std::endl;

  float boxVolume =1, foamVolume = 0;
  grid.rMax[0] = grid.rMax[0];
  for(int i=0; i<_DIMENSIONS; i++){
    boxVolume *= (grid.rMax[i]-grid.rMin[i]);
  }
  if(_DIMENSIONS == 2){
    for (int p=0; p < foam.size(); p++){
      foamVolume += M_PI*foam[p].radius*foam[p].radius;
    }
    //    foamVolume = foam.size()*M_PI*grid.radius*grid.radius;
  }
  if(_DIMENSIONS == 3){
    for (int p=0; p < foam.size(); p++){
      foamVolume += 4.0/3.0*M_PI*foam[p].radius*foam[p].radius*foam[p].radius;
    }
    //foamVolume = foam.size()*4.0/3.0*M_PI*grid.radius*grid.radius*grid.radius;
  }

  float fillingFactor = foamVolume/boxVolume;
  std::cout << "   filling factor = " << foamVolume/boxVolume << std::endl ;
  std::ofstream of1;
  of1.open("pippo.xyz", std::ofstream::out | std::ofstream::trunc);
  of1 << foam.size() << std::endl;
  of1 << "La schiuma dello Sgatto" << std::endl;
  for (int p=0; p < foam.size(); p++){
    of1 << "A ";
    for(int i=0; i<_DIMENSIONS; i++){
      of1 << foam[p].coord[i] << " ";
    }
    of1 <<  std::endl;
  }
  of1.close();

  int Npart=foam.size();
  int pointerSize = Npart*(3+1);
  float *spheresCoords = new float[pointerSize];
  for (int p=0; p < Npart; p++){
    for(int i=0; i<3; i++){
      spheresCoords[p*4+i] = foam[p].coord[i];
    }
    spheresCoords[p*4+3] = foam[p].radius;
  }
  of1.open("spheres.bin", std::ofstream::out | std::ofstream::trunc);
  of1.write((char*)&Npart,sizeof(int));
  of1.write((char*)&fillingFactor,sizeof(float));
  of1.write((char*)grid.rMin,sizeof(float)*3);
  of1.write((char*)grid.rMax,sizeof(float)*3);
  of1.write((char*)spheresCoords,sizeof(float)*pointerSize);
  of1.close();



  of1.open("spheres.txt", std::ofstream::out | std::ofstream::trunc);
  of1 << "#  " << Npart << std::endl;
  of1 << "#  " << fillingFactor << std::endl;
  of1 << "#  " << grid.rMin[0] <<  " " << grid.rMin[1] <<  " " << grid.rMin[2] <<  " " << std::endl;
  of1 << "#  " << grid.rMax[0] <<  " " << grid.rMax[1] <<  " " << grid.rMax[2] <<  " " << std::endl;
  for (int p=0; p < foam.size(); p++){
    for(int i=0; i<4; i++){
      of1 << spheresCoords[p*4+i] << " ";
    }
    of1 <<  std::endl;
  }
  of1.close();


  return 0;
}
