
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
  int k, n, i, FLAG_number = 0, FLAG_log = 0, count = 0, FLAG_angle = 0;
  int FLAG_weight=0;
  double gamma, E, m, angle, sel_angle = 180;
    double *bin, Emin, Emax, dE, norm, weight = 1, myweight;
	Dfloat ptr[NUM_COMPONENTS];
	FILE *f = fopen(args[1], "r");
	FILE *h;
	char nome[200];
	int Nbin, FLAG_P;

	if (narg < 6)
	{
		printf("-----usage:  Wenergy_Spectrum   file_name   Emin(MeV)   Emax(MeV)   Nbin   electron/protons(0/1)\n");
        printf("     options:  -log (log scale)   -num (no normalization dE) -w $WEIGHT (use WEIGHT as weight for all particles)\n");
        printf("            :  -angle $ANGLE (select particle emitted within ANGLE degrees from positive x axis\n");
        exit(0);
	}

	sscanf(args[2], "%lf", &Emin);
	sscanf(args[3], "%lf", &Emax);
	sscanf(args[4], "%i", &Nbin);
	sscanf(args[5], "%i", &FLAG_P);
	for (i = 6; i < narg; i++)
    {
        if (!strncmp(args[i], "-log", 4))
            FLAG_log = 1;
        if (!strncmp(args[i], "-num", 4)){
            FLAG_number = 1;
        }
        if (!strncmp(args[i], "-w", 2)){
            weight = atof(args[i + 1]);
            FLAG_weight=1;
        }

        if (!strncmp(args[i], "-angle", 6))
            sel_angle = atof(args[i + 1]);


    }
    sprintf(nome, "energySpectrum_%s", args[1]);
    h = fopen(nome, "w");

	bin = (double*)malloc(Nbin*sizeof(double));
	memset((void*)bin, 0, Nbin*sizeof(double));

	dE = (Emax - Emin) / Nbin;
	n = 0;
	if (FLAG_P == 1)
		m = m_proton;
	if (FLAG_P == 0)
		m = m_electron;
	while (1)
	{
		fread(ptr, sizeof(Dfloat), NUM_COMPONENTS, f);
		//fscanf(f, "%lf %lf", &angle, &gamma);
		if (feof(f)) break;
		gamma = sqrt(1 + ptr[3] * ptr[3] + ptr[4] * ptr[4] + ptr[5] * ptr[5]) - 1;
        angle = atan2(ptr[4], ptr[3]) * 180 / M_PI;
        myweight = ptr[6];
        //angle = atan2(sqrt(ptr[4]*ptr[4] + ptr[5]*ptr[5]), ptr[3]) * 180 / M_PI;
        count++;
		E = gamma*m*speed_of_light*speed_of_light*J2MeV;
        if(FLAG_weight)
            myweight=weight;
		if (E >= Emin && E < Emax && fabs(angle) < sel_angle)
		{
			k = (int) ((E - Emin) / dE);
            bin[k] += myweight;
			n++;
		}
		if (!(count % 10000))
			printf("count=%i\r", count),
			fflush(stdout);
	}
	printf("\n");

	norm = 0;
	for (k = 0; k < Nbin; k++)
		norm += bin[k] * dE;


	if (FLAG_number == 1)
		norm = 1;
	for (k = 0; k < Nbin; k++)
		bin[k] /= norm;



	for (k = 0; k < Nbin; k++)
	{
		E = Emin + dE*(k + 0.5);
		if (FLAG_log == 1)
			fprintf(h, "%e %e \n", E, log10(bin[k] + 1e-6));
		else
			fprintf(h, "%e %e \n", E, bin[k]);

		//if (bin[k]>0) fprintf(h, "%e %e \n", E, log10(bin[k]/n));
	}


}

