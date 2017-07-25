
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

#include "Wheader.h"


int main(int narg, char **args)
{
	int k, n, i, FLAG_number = 0, FLAG_log = 0, count = 0;
	double gamma, gammamin, E, m, angle;
  double *bin, thetaMin, thetaMax, dTheta, norm, weight = 1;
	Dfloat ptr[COMPONENTI];
	FILE *f = fopen(args[1], "r");
	FILE *h;
	char nome[200];
	int Nbin;

	if (narg < 6)
	{
		printf("-----usage:  angular_Spectrum   file_name   anglemin(grad)   anglemax(grad)   Nbin   gammamin\n");
		exit(0);
	}

  sscanf(args[2], "%lf", &thetaMin);
  sscanf(args[3], "%lf", &thetaMax);
	sscanf(args[4], "%i", &Nbin);
	sscanf(args[5], "%lf", &gammamin);
	for (i = 6; i < narg; i++)
	{
		if (!strncmp(args[i], "-log", 4))
			FLAG_log = 1;
		if (!strncmp(args[i], "-num", 4))
			FLAG_number = 1;
		if (!strncmp(args[i], "-w", 2))
			weight = atof(args[i + 1]);
	}
	sprintf(nome, "angularSpectrum_%s", args[1]);
	h = fopen(nome, "w");

	bin = (double*)malloc(Nbin*sizeof(double));
	memset((void*)bin, 0, Nbin*sizeof(double));

  dTheta = (thetaMax - thetaMin) / Nbin;
	n = 0;

	while (1)
	{
		fread(ptr, sizeof(Dfloat), COMPONENTI, f);
		//fscanf(f, "%lf %lf", &angle, &gamma);
		if (feof(f)) break;
		gamma = sqrt(1 + ptr[3] * ptr[3] + ptr[4] * ptr[4] + ptr[5] * ptr[5]) - 1;
    angle = atan2(ptr[4], ptr[3]) * 180 / M_PI;
		count++;

    if (angle >= thetaMin && angle<thetaMax&& gamma>gammamin)
		{
      k = (int)((angle - thetaMin) / dTheta);
      bin[k] += weight / dTheta;
			n++;
		}
		if (!(count % 10000))
			printf("count=%i\r", count),
			fflush(stdout);
	}
	printf("\n");

	norm = 0;
	for (k = 0; k < Nbin; k++)
    norm += bin[k] * dTheta;


	if (FLAG_number == 1)
		norm = 1;
	for (k = 0; k < Nbin; k++)
		bin[k] /= norm;



	for (k = 0; k < Nbin; k++)
	{
    E = thetaMin + dTheta*(k + 0.5);
		if (FLAG_log == 1)
			fprintf(h, "%e %e \n", E, log10(bin[k] + 1e-6));
		else
			fprintf(h, "%e %e \n", E, bin[k]);

		//if (bin[k]>0) fprintf(h, "%e %e \n", E, log10(bin[k]/n));
	}


}

