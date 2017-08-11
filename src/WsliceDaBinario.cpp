
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


void swap_endian_f(float *, int);

int controllo(Dfloat *, Dfloat, Dfloat);

char check_domain(Dfloat *, Dfloat [6][2]);


int main(int narg, char **args)
{
	Dfloat ptr[DIMENSIONI], *particelle, angle = 180, theta, estremi[4], EXTREM[4];
	double DDEXTREM[4];
	int i, sample = 1, np = 0, jump = 1, flag_XMIN, flag_XMAX, flag_YMIN, flag_YMAX;
	FILE *f = fopen(args[1], "r");
	FILE *g1;

	char nome[100], etichetta[100];
	int modalita, FLAG_angle = 0, FLAG_jump = 0;
	Dfloat XMIN, XMAX, YMIN, YMAX, DX, DY, mx, my;
	Dfloat **MATRIX, **MATRIX2;
	int NX = 600, NY = 600, ii, jj, n;
	int Xcomp = -1, Ycomp = -1, Zcomp = -1, Tcomp = -1; // componenti di cui si fa la proiezione
	Dfloat EXTREMS[6][2];
	size_t N, QUANTI;
	//	int index;
	int nptot, LIMITE = 10000000, NP_TOT, N_LETTURE, NP_RESTO, il;
	Dfloat factor_w = -2;
	double lettura;
	for (i = 0; i < 6; i++)
		EXTREMS[i][0] = (Dfloat)-1e30,
		EXTREMS[i][1] = (Dfloat) 1e30;
	/*
	EXTREMS[0][0]=0,
	EXTREMS[0][1]=10;
	EXTREMS[2][0]=-10,
	EXTREMS[2][1]=10;
	*/
	flag_XMIN = flag_XMAX = flag_YMIN = flag_YMAX = 0;
	if (narg <= 1)
	{
		printf("\nusage: WsliceDaBinario  file_in  $MODALITY\n");
		printf("\nMODALITY:\n");
		printf("\t0=XY, 1=XPX, 2=ZPZ, 3=YPX, 4=theta_gamma-1,\n");
		printf("\t5=theta_MeV, 6=ZY_gamma-1, 7=ZX_PX, 8=XY_gamma-1, 9=PZPX\n");
		printf("\t10=PYPZ, 11=ZX_JX, 12=ZX_JZ, 13=XY_MeV,  14=ZY\n");
		printf("\t15=PYPZ_PX, 16=XZ, 17=YPY, 18=Ytheta\n");
		printf("OPTIONS:\n");
		printf("\t-sel $DEGREES\n");
		printf("\t-jump $JUMPS (sub sampling)\n");
		printf("\t-nx $NXBIN -ny $NYBIN\n");
		printf("\t-w $FORCED_WEIGHT\n");
		printf("   binning extrema: \n");
		printf("\t$XMIN^$XMAX ($XMIN^ or ^$XMAX are also allowed)\n");
		printf("\t$YMIN#$YMAX ($YMIN# or #$YMAX are also allowed)\n");
		printf("   extrema of the selection:\n");
		printf("\t-xmin $XMIN -xmax $XMAX -ymin $YMIN -ymax $YMAX -zmin $ZMIN -zmax $ZMAX\n");
		printf("\t-pxmin $PXMIN -pxmax $PXMAX -pymin $PYMIN -pymax $PYMAX -pzmin $PZMIN -pzmax $PZMAX\n\n");
		return 123;
	}
	else
	{
		sscanf(args[2], "%i", &modalita);
		for (i = 3; i < narg; i++)
		{
			if (!strncmp(args[i], "-sel", 4))
			{
				FLAG_angle = 1;
				sscanf(args[i + 1], "%lf", &lettura);
				angle = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-nx", 3))
			{
				sscanf(args[i + 1], "%i", &NX);
			}
			if (!strncmp(args[i], "-ny", 3))
			{
				sscanf(args[i + 1], "%i", &NY);
			}
			if (!strncmp(args[i], "-zmin", 4))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[2][0] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-zmax", 5))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[2][1] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-ymin", 4))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[1][0] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-ymax", 5))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[1][1] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-xmin", 5))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[0][0] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-xmax", 5))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[0][1] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-pxmin", 6))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[3][0] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-pxmax", 6))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[3][1] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-pymin", 6))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[4][0] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-pymax", 6))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[4][1] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-pzmin", 6))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[5][0] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-pzmax", 6))
			{
				sscanf(args[i + 1], "%lf", &lettura);
				EXTREMS[5][1] = (Dfloat)lettura;
			}
			if (!strncmp(args[i], "-jump", 5))
			{
				FLAG_jump = 1;
				sscanf(args[i + 1], "%i", &jump);
			}
			if (args[i][0] == '^')
			{
				flag_XMIN = FALSE;
				flag_XMAX = TRUE;
				sscanf(args[i], "^%lf", &DDEXTREM[2]);
				EXTREM[2] = (Dfloat)DDEXTREM[2];
			}
			else if (args[i][strlen(args[i]) - 1] == '^')
			{
				flag_XMIN = TRUE;
				flag_XMAX = FALSE;
				sscanf(args[narg - i], "%lf^", &DDEXTREM[0]);
				EXTREM[0] = (Dfloat)DDEXTREM[0];
			}
			else if (strstr(args[i], "^"))
			{
				flag_XMIN = TRUE;
				flag_XMAX = TRUE;
				sscanf(args[i], "%lf^%lf", &DDEXTREM[0], &DDEXTREM[2]);
				EXTREM[2] = (Dfloat)DDEXTREM[2];
				EXTREM[0] = (Dfloat)DDEXTREM[0];
			}


			//ora estremi per la y
			if (args[i][0] == '#')
			{
				flag_YMIN = FALSE;
				flag_YMAX = TRUE;
				sscanf(args[i], "#%lf", &DDEXTREM[3]);
				EXTREM[3] = (Dfloat)DDEXTREM[3];
			}
			else if (args[i][strlen(args[i]) - 1] == '#')
			{
				flag_YMIN = TRUE;
				flag_YMAX = FALSE;
				sscanf(args[i], "%lf#", &DDEXTREM[1]);
				EXTREM[1] = (Dfloat)DDEXTREM[1];
			}
			else if (strstr(args[i], "#"))
			{
				flag_YMIN = TRUE;
				flag_YMAX = TRUE;
				sscanf(args[i], "%lf#%lf", &DDEXTREM[1], &DDEXTREM[3]);
				EXTREM[1] = (Dfloat)DDEXTREM[1];
				EXTREM[3] = (Dfloat)DDEXTREM[3];
			}
		}
	}
	for (i = 1; i < narg; i++)
		if (!strncmp(args[i], "-w", 2))
			factor_w = (Dfloat) (atof(args[i + 1]));

	estremi[0] = estremi[1] = (Dfloat)-1e+37;
	estremi[2] = estremi[3] = (Dfloat)+1e+37;

	printf("FLAG_jump=%i  jump=%i\n", FLAG_jump, jump);
	printf("FLAG_angle=%i  angle=%f\n", FLAG_angle, angle);
	printf("estremi  %g  %g  %g  %g\n", estremi[0], estremi[1], estremi[2], estremi[3]);
	fflush(stdout);

	switch (modalita)
	{
	case 0:
		sprintf(nome, "SLICE_%s_XY", args[1]);
		sprintf(etichetta, "density_XY");
		g1 = fopen(nome, "w");
		Xcomp = 0;
		Ycomp = 1;
		break;
	case 1:
		sprintf(nome, "SLICE_%s_XPX", args[1]);
		sprintf(etichetta, "density_XPX");
		g1 = fopen(nome, "w");
		Xcomp = 0;
		Ycomp = 3;
		break;
	case 2:
		sprintf(nome, "SLICE_%s_ZPZ", args[1]);
		sprintf(etichetta, "density_ZPZ");
		g1 = fopen(nome, "w");
		Xcomp = 2;
		Ycomp = 5;
		break;
	case 3:
		sprintf(nome, "SLICE_%s_YPX", args[1]);
		sprintf(etichetta, "density_YPZ");
		g1 = fopen(nome, "w");
		Xcomp = 1;
		Ycomp = 3;
		break;
	case 4:
		sprintf(nome, "SLICE_%s_theta_gamma-1", args[1]);
		sprintf(etichetta, "density_ZPZ");
		g1 = fopen(nome, "w");
		Xcomp = THETA;
		Ycomp = GAMMA1;
		break;
	case 5:
		sprintf(nome, "SLICE_%s_theta_MeV", args[1]);
		g1 = fopen(nome, "w");
		Xcomp = THETA;
		Ycomp = ENERGIA;
		break;
	case 6:
		sprintf(nome, "SLICE_%s_ZY_gamma-1", args[1]);
		g1 = fopen(nome, "w");
		Xcomp = 2;
		Ycomp = 1;
		Zcomp = GAMMA1;
		break;
	case 7:
		sprintf(nome, "SLICE_%s_ZX_PX", args[1]);
		sprintf(etichetta, "density_ZPZ");
		g1 = fopen(nome, "w");
		Xcomp = 2;
		Ycomp = 0;
		Zcomp = 3;
		printf("Zcomp=%i\n", Zcomp);
		break;
	case 8:
		sprintf(nome, "SLICE_%s_XY_gamma-1", args[1]);
		sprintf(etichetta, "density_ZPZ");
		g1 = fopen(nome, "w");
		Xcomp = 0;
		Ycomp = 1;
		Zcomp = GAMMA1;
		break;
	case 9:
		sprintf(nome, "SLICE_%s_PZPX", args[1]);
		g1 = fopen(nome, "w");
		Xcomp = 5;
		Ycomp = 3;
		break;
	case 10:
		sprintf(nome, "SLICE_%s_PYPZ", args[1]);
		g1 = fopen(nome, "w");
		Xcomp = 4;
		Ycomp = 5;
		break;
	case 11:
		sprintf(nome, "SLICE_%s_ZX_JX", args[1]);
		sprintf(etichetta, "density_ZPZ");
		g1 = fopen(nome, "w");
		Xcomp = 2;
		Ycomp = 0;
		Tcomp = VX;
		break;
	case 12:
		sprintf(nome, "SLICE_%s_ZX_JZ", args[1]);
		sprintf(etichetta, "density_ZPZ");
		g1 = fopen(nome, "w");
		Xcomp = 2;
		Ycomp = 0;
		Tcomp = VZ;
		break;
	case 13:
		sprintf(nome, "SLICE_%s_XY_MeV", args[1]);
		sprintf(etichetta, "density_ZPZ");
		g1 = fopen(nome, "w");
		Xcomp = 0;
		Ycomp = 1;
		Zcomp = ENERGIA;
		break;
	case 14:

		sprintf(nome, "SLICE_%s_ZY", args[1]);
		sprintf(etichetta, "density_ZPZ");
		g1 = fopen(nome, "w");
		Xcomp = 2;
		Ycomp = 1;
		break;
	case 15:
		sprintf(nome, "SLICE_%s_PYPZ_PX", args[1]);
		sprintf(etichetta, "density_ZPZ");
		g1 = fopen(nome, "w");
		Xcomp = 4;
		Ycomp = 5;
		Zcomp = 3;
		break;
	case 16:
		sprintf(nome, "SLICE_%s_XZ", args[1]);
		sprintf(etichetta, "density_XZ");
		g1 = fopen(nome, "w");
		Xcomp = 0;
		Ycomp = 2;
		break;
	case 17:
		sprintf(nome, "SLICE_%s_YPY_PX", args[1]);
		sprintf(etichetta, "density_YPY");
		g1 = fopen(nome, "w");
		Xcomp = 1;
		Ycomp = 4;
		break;
	case 18:
		sprintf(nome, "SLICE_%s_Ytheta", args[1]);
		sprintf(etichetta, "density_Ytheta");
		g1 = fopen(nome, "w");
		Xcomp = 1;
		Ycomp = THETA;
		break;
	default:
		printf("no valid modality has been given");
		return 117;
	}


	printf("Inizio letture \n");
	fflush(stdout);
	fseeko(f, 0, SEEK_END);
	N = (size_t)ftello(f);
	rewind(f);
	QUANTI = (N / sizeof(Dfloat));
	NP_TOT = (int)(QUANTI / NUM_COMPONENTS);
	printf("Il file contiene %i particelle \n", NP_TOT);
	fflush(stdout);

	N_LETTURE = NP_TOT / LIMITE;
	NP_RESTO = NP_TOT%LIMITE;
	nptot = LIMITE;   //HO SETTATO nptot=LIMITE!!!!!
	printf("occorre effettuare %i letture di %i particelle e\n", N_LETTURE, LIMITE);
	printf("\t una lettura di %i particelle\n", NP_RESTO);
	fflush(stdout);

	particelle = (Dfloat*)malloc(sizeof(Dfloat)*LIMITE*NUM_COMPONENTS);
	QUANTI = LIMITE*NUM_COMPONENTS;

	XMIN = (Dfloat) 1e20;
	XMAX = (Dfloat)-1e20;
	YMIN = (Dfloat) 1e20;
	YMAX = (Dfloat)-1e20;
	printf("PRIMA:\nXMIN=%g\tXMAX=%g\tYMIN=%g\tYMAX=%g\n", XMIN, XMAX, YMIN, YMAX);
	printf("Xcomp = %i    Ycomp = %i\n", Xcomp, Ycomp);
	if (!(flag_XMIN&&flag_YMIN&&flag_XMAX&&flag_YMAX))
	{
		printf("devo leggere per trovare estremi\n");
		printf("flag_XMIN=%i\tflag_YMIN=%i\tflag_XMAX=%i\tflag_YMAX=%i\n", flag_XMIN, flag_YMIN, flag_XMAX, flag_YMAX);
		for (il = 0; il < N_LETTURE; il++)
		{
			printf("inizio lettura numero %i\n", il);
			fflush(stdout);
			nptot = LIMITE;   //HO SETTATO nptot=LIMITE!!!!!
			QUANTI = nptot*NUM_COMPONENTS;
			fread(particelle, sizeof(Dfloat), QUANTI, f);
			for (n = 0; n < nptot; n++)
			{
				for (i = 0; i < NUM_COMPONENTS; i++) ptr[i] = particelle[n*NUM_COMPONENTS + i];
				ptr[GAMMA1] = (Dfloat)(sqrt(1. + ptr[3] * ptr[3] + ptr[4] * ptr[4] + ptr[5] * ptr[5]) - 1.);
				theta = (Dfloat)(atan2(sqrt(ptr[4] * ptr[4] + ptr[5] * ptr[5]), ptr[3]) * 180. / M_PI);
				ptr[THETA] = theta;
				ptr[ENERGIA] = (Dfloat)(ptr[GAMMA1] * m_proton_mev);
				ptr[VX] = (Dfloat)(ptr[3] / (ptr[GAMMA1] + 1.));
				ptr[VZ] = (Dfloat)(ptr[5] / (ptr[GAMMA1] + 1.));

				mx = ptr[Xcomp];
				my = ptr[Ycomp];

				XMIN = MIN(XMIN, mx);
				YMIN = MIN(YMIN, my);
				XMAX = MAX(XMAX, mx);
				YMAX = MAX(YMAX, my);
			}
		}
		printf("inizio lettura del RESTO\n");
		fflush(stdout);

		nptot = NP_RESTO;
		QUANTI = nptot*NUM_COMPONENTS;
		fread(particelle, sizeof(Dfloat), QUANTI, f);
		for (n = 0; n < nptot; n++)
		{
			for (i = 0; i < NUM_COMPONENTS; i++) ptr[i] = particelle[n*NUM_COMPONENTS + i];
			ptr[GAMMA1] = (Dfloat)(sqrt(1. + ptr[3] * ptr[3] + ptr[4] * ptr[4] + ptr[5] * ptr[5]) - 1.);
			theta = (Dfloat)(atan2(sqrt(ptr[4] * ptr[4] + ptr[5] * ptr[5]), ptr[3]) * 180. / M_PI);
			ptr[THETA] = theta;
			ptr[ENERGIA] = (Dfloat)(ptr[GAMMA1] * m_proton_mev);
			ptr[VX] = (Dfloat)(ptr[3] / (ptr[GAMMA1] + 1.));
			ptr[VZ] = (Dfloat)(ptr[5] / (ptr[GAMMA1] + 1.));

			mx = ptr[Xcomp];
			my = ptr[Ycomp];

			XMIN = MIN(XMIN, mx);
			YMIN = MIN(YMIN, my);
			XMAX = MAX(XMAX, mx);
			YMAX = MAX(YMAX, my);
		}
		// ---------------------------
		//   RIAVVOLGO LO STREAM!!!!
		// ---------------------------
		rewind(f);
	}
	// ---------------------------

	if (flag_XMIN)	XMIN = EXTREM[0];
	if (flag_YMIN)	YMIN = EXTREM[1];
	if (flag_XMAX)	XMAX = EXTREM[2];
	if (flag_YMAX)	YMAX = EXTREM[3];

	printf("NX=%i  NY=%i\n", NX, NY);
	printf("DOPO:\nXMIN=%g\tXMAX=%g\tYMIN=%g\tYMAX=%g\n", XMIN, XMAX, YMIN, YMAX);

	MATRIX = (Dfloat**)malloc(NX*sizeof(Dfloat*));
	for (ii = 0; ii < NX; ii++)
	{
		MATRIX[ii] = (Dfloat*)malloc(NY*sizeof(Dfloat));
		for (jj = 0; jj < NY; jj++) MATRIX[ii][jj] = 0;
	}
	DX = (XMAX - XMIN) / (Dfloat)NX;
	DY = (YMAX - YMIN) / (Dfloat)NY;
	printf("DX=%f  DY=%f\n", DX, DY);
	printf("numero particelle: %i\n", nptot);
	printf("Xcomp = %i\tYcomp = %i\tZcomp = %i\tTcomp = %i\n", Xcomp, Ycomp, Zcomp, Tcomp);

	if (Zcomp < 0 && Tcomp < 0)
	{
		for (il = 0; il < N_LETTURE; il++)
		{
			printf("inizio lettura numero %i\n", il);
			fflush(stdout);

			nptot = LIMITE;   //HO SETTATO nptot=LIMITE!!!!!
			QUANTI = nptot*NUM_COMPONENTS;
			fread(particelle, sizeof(Dfloat), QUANTI, f);

			printf("ho finito la lettura numero %i, inizio calcoli\n", il);
			fflush(stdout);
			for (n = 0; n < nptot; n++)
			{
				for (i = 0; i < NUM_COMPONENTS; i++) ptr[i] = particelle[n*NUM_COMPONENTS + i];
				ptr[GAMMA1] = (Dfloat)(sqrt(1. + ptr[3] * ptr[3] + ptr[4] * ptr[4] + ptr[5] * ptr[5]) - 1.);
				theta = (Dfloat)(atan2(sqrt(ptr[4] * ptr[4] + ptr[5] * ptr[5]), ptr[3]) * 180. / M_PI);
				ptr[THETA] = theta;
				ptr[ENERGIA] = (Dfloat)(ptr[GAMMA1] * m_proton_mev);
				ptr[VX] = (Dfloat)(ptr[3] / (ptr[GAMMA1] + 1.));
				ptr[VZ] = (Dfloat)(ptr[5] / (ptr[GAMMA1] + 1.));

				mx = ptr[Xcomp];
				my = ptr[Ycomp];
				if (mx >= XMIN&&mx <= XMAX&&my >= YMIN&&my <= YMAX)
				{
					ii = (int)MIN((mx - XMIN) / DX, NX - 1);
					jj = (int)MIN((my - YMIN) / DY, NY - 1);
					//printf("ii=%i jj=%i\n",ii, jj);
					//fflush(stdout);
					if (check_domain(ptr, EXTREMS))
					{
						if (factor_w < 0)	MATRIX[ii][jj] += (Dfloat)(fabs(ptr[6])*1e8 / (DX*DY));
						else			MATRIX[ii][jj] += (Dfloat)(factor_w);
					}
				}
			}
		}

		printf("inizio lettura del RESTO\n");
		fflush(stdout);
		nptot = NP_RESTO;
		QUANTI = nptot*NUM_COMPONENTS;
		fread(particelle, sizeof(Dfloat), QUANTI, f);

		for (n = 0; n < nptot; n++)
		{
			for (i = 0; i < NUM_COMPONENTS; i++) ptr[i] = particelle[n*NUM_COMPONENTS + i];
			ptr[GAMMA1] = (Dfloat)(sqrt(1. + ptr[3] * ptr[3] + ptr[4] * ptr[4] + ptr[5] * ptr[5]) - 1.);
			theta = (Dfloat)(atan2(sqrt(ptr[4] * ptr[4] + ptr[5] * ptr[5]), ptr[3]) * 180. / M_PI);
			ptr[THETA] = theta;
			ptr[ENERGIA] = (Dfloat)(ptr[GAMMA1] * m_proton_mev);
			ptr[VX] = (Dfloat)(ptr[3] / (ptr[GAMMA1] + 1.));
			ptr[VZ] = (Dfloat)(ptr[5] / (ptr[GAMMA1] + 1.));

			mx = ptr[Xcomp];
			my = ptr[Ycomp];
			if (mx >= XMIN&&mx <= XMAX&&my >= YMIN&&my <= YMAX)
			{
				ii = (int)MIN((mx - XMIN) / DX, NX - 1);
				jj = (int)MIN((my - YMIN) / DY, NY - 1);
				if (check_domain(ptr, EXTREMS))
				{
					if (factor_w<0)	MATRIX[ii][jj] += (Dfloat)(fabs(ptr[6])*1e8 / (DX*DY));
					else			MATRIX[ii][jj] += (Dfloat)(factor_w);
				}
			}
		}
	}
	else if (Zcomp>0 && Tcomp < 0)
	{
		MATRIX2 = (Dfloat**)malloc(NX*sizeof(Dfloat*));
		for (ii = 0; ii < NX; ii++)
		{
			MATRIX2[ii] = (Dfloat*)malloc(NY*sizeof(Dfloat));
			for (jj = 0; jj < NY; jj++)
				MATRIX2[ii][jj] = 0;
		}
		for (il = 0; il < N_LETTURE; il++)
		{
			printf("inizio lettura numero %i\n", il);
			fflush(stdout);
			nptot = LIMITE;   //HO SETTATO nptot=LIMITE!!!!!
			QUANTI = nptot*NUM_COMPONENTS;
			fread(particelle, sizeof(Dfloat), QUANTI, f);
			printf("ho finito la lettura numero %i, inizio calcoli\n", il);
			fflush(stdout);

			for (n = 0; n < nptot; n++)
			{
				for (i = 0; i < NUM_COMPONENTS; i++) ptr[i] = particelle[n*NUM_COMPONENTS + i];
				ptr[GAMMA1] = (Dfloat)(sqrt(1. + ptr[3] * ptr[3] + ptr[4] * ptr[4] + ptr[5] * ptr[5]) - 1.);
				theta = (Dfloat)(atan2(sqrt(ptr[4] * ptr[4] + ptr[5] * ptr[5]), ptr[3]) * 180. / M_PI);
				ptr[THETA] = theta;
				ptr[ENERGIA] = (Dfloat)(ptr[GAMMA1] * m_proton_mev);
				ptr[VX] = (Dfloat)(ptr[3] / (ptr[GAMMA1] + 1.));
				ptr[VZ] = (Dfloat)(ptr[5] / (ptr[GAMMA1] + 1.));

				mx = ptr[Xcomp];
				my = ptr[Ycomp];
				if (mx >= XMIN&&mx < XMAX&&my >= YMIN&&my < YMAX)
				{
					ii = (int)((mx - XMIN) / DX);
					jj = (int)((my - YMIN) / DY);
					if (check_domain(ptr, EXTREMS))
					{
						MATRIX[ii][jj] += ptr[6] * ptr[Zcomp];
						MATRIX2[ii][jj] += ptr[6];
					}
				}
			}

		}

		printf("inizio lettura del RESTO\n");
		fflush(stdout);
		nptot = NP_RESTO;
		QUANTI = nptot*NUM_COMPONENTS;
		fread(particelle, sizeof(Dfloat), QUANTI, f);

		for (n = 0; n < nptot; n++)
		{
			for (i = 0; i < NUM_COMPONENTS; i++) ptr[i] = particelle[n*NUM_COMPONENTS + i];
			ptr[GAMMA1] = (Dfloat)(sqrt(1. + ptr[3] * ptr[3] + ptr[4] * ptr[4] + ptr[5] * ptr[5]) - 1.);
			theta = (Dfloat)(atan2(sqrt(ptr[4] * ptr[4] + ptr[5] * ptr[5]), ptr[3]) * 180. / M_PI);
			ptr[THETA] = theta;
			ptr[ENERGIA] = (Dfloat)(ptr[GAMMA1] * m_proton_mev);
			ptr[VX] = (Dfloat)(ptr[3] / (ptr[GAMMA1] + 1.));
			ptr[VZ] = (Dfloat)(ptr[5] / (ptr[GAMMA1] + 1.));

			mx = ptr[Xcomp];
			my = ptr[Ycomp];
			if (mx >= XMIN&&mx < XMAX&&my >= YMIN&&my < YMAX)
			{
				ii = (int)((mx - XMIN) / DX);
				jj = (int)((my - YMIN) / DY);
				if (check_domain(ptr, EXTREMS))
				{
					MATRIX[ii][jj] += ptr[6] * ptr[Zcomp];
					MATRIX2[ii][jj] += ptr[6];
				}
			}
		}
		for (ii = 0; ii < NX; ii++)
			for (jj = 0; jj < NY; jj++)
			{
				//printf("i=%i j=%i  num=%e denum=%e\n", ii, jj,MATRIX[ii][jj],MATRIX2[ii][jj]  );
				if (MATRIX2[ii][jj] != 0)
					MATRIX[ii][jj] /= MATRIX2[ii][jj];
			}
	}
	else
	{
		for (il = 0; il < N_LETTURE; il++)
		{
			printf("inizio lettura numero %i\n", il);
			fflush(stdout);

			nptot = LIMITE;   //HO SETTATO nptot=LIMITE!!!!!
			QUANTI = nptot*NUM_COMPONENTS;
			fread(particelle, sizeof(Dfloat), QUANTI, f);

			printf("ho finito la lettura numero %i, inizio calcoli\n", il);
			fflush(stdout);
			for (n = 0; n < nptot; n++)
			{
				for (i = 0; i < NUM_COMPONENTS; i++) ptr[i] = particelle[n*NUM_COMPONENTS + i];
				ptr[GAMMA1] = (Dfloat)(sqrt(1. + ptr[3] * ptr[3] + ptr[4] * ptr[4] + ptr[5] * ptr[5]) - 1.);
				theta = (Dfloat)(atan2(sqrt(ptr[4] * ptr[4] + ptr[5] * ptr[5]), ptr[3]) * 180. / M_PI);
				ptr[THETA] = theta;
				ptr[ENERGIA] = (Dfloat)(ptr[GAMMA1] * m_proton_mev);
				ptr[VX] = (Dfloat)(ptr[3] / (ptr[GAMMA1] + 1.));
				ptr[VZ] = (Dfloat)(ptr[5] / (ptr[GAMMA1] + 1.));

				mx = ptr[Xcomp];
				my = ptr[Ycomp];
				if (mx >= XMIN&&mx <= XMAX&&my >= YMIN&&my <= YMAX)
				{
					ii = (int)MIN((mx - XMIN) / DX, NX - 1);
					jj = (int)MIN((my - YMIN) / DY, NY - 1);
					//printf("ii=%i jj=%i\n",ii, jj);
					//fflush(stdout);
					if (check_domain(ptr, EXTREMS))
					{
						if (factor_w < 0)	MATRIX[ii][jj] += (Dfloat)(ptr[6] * ptr[Tcomp] * 1e8 / (DX*DY));
						else			MATRIX[ii][jj] += (Dfloat)(factor_w*ptr[Tcomp]);
					}
				}
			}
		}

		printf("inizio lettura del RESTO\n");
		fflush(stdout);
		nptot = NP_RESTO;
		QUANTI = nptot*NUM_COMPONENTS;
		fread(particelle, sizeof(Dfloat), QUANTI, f);

		for (n = 0; n < nptot; n++)
		{
			for (i = 0; i < NUM_COMPONENTS; i++) ptr[i] = particelle[n*NUM_COMPONENTS + i];
			ptr[GAMMA1] = (Dfloat)(sqrt(1. + ptr[3] * ptr[3] + ptr[4] * ptr[4] + ptr[5] * ptr[5]) - 1.);
			theta = (Dfloat)(atan2(sqrt(ptr[4] * ptr[4] + ptr[5] * ptr[5]), ptr[3]) * 180. / M_PI);
			ptr[THETA] = theta;
			ptr[ENERGIA] = (Dfloat)(ptr[GAMMA1] * m_proton_mev);
			ptr[VX] = (Dfloat)(ptr[3] / (ptr[GAMMA1] + 1.));
			ptr[VZ] = (Dfloat)(ptr[5] / (ptr[GAMMA1] + 1.));

			mx = ptr[Xcomp];
			my = ptr[Ycomp];
			if (mx >= XMIN&&mx <= XMAX&&my >= YMIN&&my <= YMAX)
			{
				ii = (int)MIN((mx - XMIN) / DX, NX - 1);
				jj = (int)MIN((my - YMIN) / DY, NY - 1);
				if (check_domain(ptr, EXTREMS))
				{
					if (factor_w < 0)	MATRIX[ii][jj] += (Dfloat)(ptr[6] * ptr[Tcomp] * 1e8 / (DX*DY));
					else			MATRIX[ii][jj] += (Dfloat)(factor_w*ptr[Tcomp]);
				}
			}
		}


	}
	fprintf(g1, "%i\n%i\n%i\n", NX, NY, 1);
	fprintf(g1, "%f  %f\n%f  %f\n", XMIN, YMIN, XMAX, YMAX);
	printf("%i\n%i\n%i\n", NX, NY, 1);
	printf("%f  %f\n%f  %f\n", XMIN, YMIN, XMAX, YMAX);

	for (ii = 0; ii < NX; ii++)
	{
		for (jj = 0; jj < NY; jj++)
		{
			mx = (Dfloat)(XMIN + DX*(ii + 0.5));
			my = (Dfloat)(YMIN + DY*(jj + 0.5));
			fprintf(g1, "%f  %f  %.6g\n", mx, my, MATRIX[ii][jj]);
		}
	}


	// VTK format
	/*
	sprintf(nome,"%s.vtk",nome);
	f1=fopen(nome, "w");

	fprintf(g1,"# vtk DataFile Version 2.0\n");
	fprintf(g1,"titolo mio\n");
	fprintf(g1,"ASCII\n");
	fprintf(g1,"DATASET STRUCTURED_POINTS\n");
	fprintf(g1,"DIMENSIONS %i %i\n",NX, NY);
	fprintf(g1,"ORIGIN %f %f\n",XMIN, YMIN);
	dx=(xmax-xmin)/(nx1-1);
	dy=(ymax-ymin)/(ny1-1);
	dz=(zmax-zmin)/(nz1-1);
	fprintf(g1,"SPACING %f %f\n",DX, DY);
	fprintf(g1,"POINT_DATA %i\n",NX*NY);
	fprintf(g1,"%s \n",etichetta);
	fprintf(g1,"LOOKUP_TABLE default\n");
	for(ii=0;ii<NX;ii++)
	{
	for(jj=0;jj<NY;jj++)
	{
	fprintf(g1, "%.6g\n",MATRIX[ii][jj]);
	}
	}
	*/

	fclose(g1);
}



void swap_endian_f(float* in_f, int n)
{
	int i;

	union
	{
		int imio;
		float fmio;
		char arr[4];
	}x;

	char buff;
	for (i = 0; i < n; i++)
	{
		x.fmio = in_f[i];

		buff = x.arr[0];
		x.arr[0] = x.arr[3];
		x.arr[3] = buff;
		buff = x.arr[1];
		x.arr[1] = x.arr[2];
		x.arr[2] = buff;

		in_f[i] = x.fmio;
	}
}


int controllo(Dfloat *estremi, Dfloat X, Dfloat Y)
{
	if (X <= estremi[2] && X >= estremi[0] && Y <= estremi[3] && Y >= estremi[1])
		return 1;
	else
		return 0;

}


char check_domain(Dfloat *X, Dfloat EXTREMS[][2])
{
	int  i;
	char IN = 1;

	for (i = 0; i < 6; i++)
		IN *= (X[i] >= EXTREMS[i][0] && X[i] < EXTREMS[i][1]);

	return IN;
}


