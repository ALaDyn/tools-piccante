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

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

#define m_electron 9.1095e-31      //electrons
#define m_proton 1.6726231e-27   //protons
#define c 299792458
#define J2MeV (1e-6/1.60218e-19)
#define COMPONENTI 7
#define Dfloat float 

//#define N 250

int main(int narg, char **args)
{
  int k, n, i,FLAG_number=0, FLAG_log=0, count=0;
  double z, uz, ux, uy, x, y, gamma, gammamin, E, m, angle;
  double *bin, Emin, Emax, dE, norm, weight=1;
  Dfloat ptr[COMPONENTI];
  FILE *f=fopen(args[1], "r");
  FILE *h;
  char nome[200];
  int Nbin, FLAG_P;
  
  if(narg<6)
    {
      printf("-----usage:  Wenergy_Spectrum   file_name   Emin(MeV)   Emax(MeV)   Nbin   electron/protons(0/1)\n");
      exit(0);
    }
    
  sscanf(args[2], "%lf", &Emin);
  sscanf(args[3], "%lf", &Emax);
  sscanf(args[4], "%i", &Nbin);
  sscanf(args[5], "%i", &FLAG_P);
  for (i=6; i<narg; i++)
     {
          if (!strncmp(args[i], "-log",4))
              FLAG_log=1;
          if (!strncmp(args[i], "-num",4))
              FLAG_number=1;
	  if (!strncmp(args[i], "-w",2))
	    weight=atof(args[i+1]);


     }
  sprintf(nome, "energySpectrum_%s", args[1]);
  h=fopen(nome, "w");

  bin= (double*)malloc(Nbin*sizeof(double));
  memset((void*)bin, 0, Nbin*sizeof(double));

  dE=(Emax-Emin)/Nbin;
  n=0;
  if(FLAG_P==1)
    m=m_proton;
  if(FLAG_P==0)
    m=m_electron;
  while(1)
    {
      fread(ptr, sizeof(Dfloat), COMPONENTI, f);
      //fscanf(f, "%lf %lf", &angle, &gamma);
      if (feof(f)) break;
      gamma=sqrt(1+ptr[3]*ptr[3]+ptr[4]*ptr[4]+ptr[5]*ptr[5])-1;
      angle=atan2(ptr[5]/gamma,ptr[3]/gamma)*180/M_PI;
      count++;
      E=gamma*m*c*c*J2MeV;

      if (E>=Emin && E<Emax)
	{
	  k=(E-Emin)/dE;
	  bin[k]+=weight/dE;
	  n++;
	}
      if(!(count%10000))
	printf("count=%i\r",count),
	  fflush(stdout);
    }
  printf("\n");

  norm=0;
  for (k=0; k<Nbin; k++)
    norm+=bin[k]*dE;
  

  if(FLAG_number==1)
    norm=1;
  for (k=0; k<Nbin; k++)
    bin[k]/=norm;



  for (k=0; k<Nbin; k++)
    {
      E=Emin+dE*(k+0.5);
      if(FLAG_log==1)
          fprintf(h, "%e %e \n", E, log10(bin[k]+1e-6));
      else
          fprintf(h, "%e %e \n", E, bin[k]);

      //if (bin[k]>0) fprintf(h, "%e %e \n", E, log10(bin[k]/n));
    }
  

}
