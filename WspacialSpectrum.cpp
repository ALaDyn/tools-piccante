
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

#include "Wheader.h"


int main(int narg, char **args)
{
  int k, n, i, FLAG_number = 0, FLAG_log = 0, count = 0, coord;
  double gamma, gammamin, E, m, angle, anglemin, anglemax;
  double *bin, rMin, rMax, deltaR, norm, weight = 1;
  Dfloat ptr[COMPONENTI];
  FILE *f = fopen(args[1], "r");
  FILE *h;
  char nome[200];
  int Nbin;

  if (narg < 6)
  {
    printf("-----usage:  spacialSpectrum   file_name   rmin   rmax coords(0,1,2)   Nbin   gammamin\n");
    printf("-----option:  -anglemin $ANGLEMIN -anglemax $ANGLEMAX\n");
    exit(0);
  }

  sscanf(args[2], "%lf", &rMin);
  sscanf(args[3], "%lf", &rMax);
  sscanf(args[4], "%i", &coord);
  sscanf(args[5], "%i", &Nbin);
  sscanf(args[6], "%lf", &gammamin);
  for (i = 6; i < narg; i++)
  {
    if (!strncmp(args[i], "-log", 4))
      FLAG_log = 1;
    if (!strncmp(args[i], "-num", 4))
      FLAG_number = 1;
    if (!strncmp(args[i], "-w", 2))
      weight = atof(args[i + 1]);
    if (!strncmp(args[i], "-anglemin", 9))
      anglemin = atof(args[i + 1]);
    if (!strncmp(args[i], "-anglemax", 9))
      anglemax = atof(args[i + 1]);
  }
  char coordName[3][2]={"x","y","z"};
  sprintf(nome, "spacialSpectrum_%s_%s", coordName[coord], args[1]);
  h = fopen(nome, "w");

  bin = (double*)malloc(Nbin*sizeof(double));
  memset((void*)bin, 0, Nbin*sizeof(double));

  deltaR = (rMax - rMin) / Nbin;
  n = 0;

  while (1)
  {
    fread(ptr, sizeof(Dfloat), COMPONENTI, f);
    //fscanf(f, "%lf %lf", &angle, &gamma);
    if (feof(f)) break;
    gamma = sqrt(1 + ptr[3] * ptr[3] + ptr[4] * ptr[4] + ptr[5] * ptr[5]) - 1;
    angle = atan2(ptr[4], ptr[3]) * 180 / M_PI;
    count++;

    if (ptr[coord] >= rMin && ptr[coord]<rMax&& gamma>gammamin)
    {
      k = (int)((angle - rMin) / deltaR);
      bin[k] += weight / deltaR;
      n++;
    }
    if (!(count % 10000))
      printf("count=%i\r", count),
      fflush(stdout);
  }
  printf("\n");

  norm = 0;
  for (k = 0; k < Nbin; k++)
    norm += bin[k] * deltaR;



  for (k = 0; k < Nbin; k++)
  {
    E = rMin + deltaR*(k + 0.5);
    if (FLAG_log == 1)
      fprintf(h, "%e %e \n", E, log10(bin[k] + 1e-6));
    else
      fprintf(h, "%e %e \n", E, bin[k]);

    //if (bin[k]>0) fprintf(h, "%e %e \n", E, log10(bin[k]/n));
  }


}

