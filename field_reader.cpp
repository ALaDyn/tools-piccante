#include<stdio.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<malloc.h>
#include<cmath>
#include<iomanip>
#include<string>
#include <stdint.h>

using namespace std;


int is_big_endian(void)
{
	union {
		uint32_t i;
		char c[4];
	} bint = {0x01020304};

	return bint.c[0] == 1; 
}

inline void drawLoadBar(long i, long Ntot, long R, int sizeBar){
	if (i % (Ntot/R) != 0) return;
	
	std::cout << "\r";

	float ratio = i /((float)Ntot);

	int numSymbols = (int) sizeBar*ratio;

	std::cout << std::setw(3) << (int)(ratio*100) << "% [";
	
	int j;
	for (j = 0; j < numSymbols; j++){
		std::cout << "=";
	}
	for (; j < sizeBar; j++){
		std::cout << " ";
	}

	std::cout <<"]";
	std::cout.flush();
}


int main (const int argc, const char *argv[])
{
	int Ncells[3], rNproc[3], locOrigin[3],locNcells[3];
	double *xiCoords,*yiCoords,*ziCoords;
	double *xhCoords,*yhCoords,*zhCoords;
	double *riCoords[3], *rhCoords[3];
	int Nproc;
	int Ncomp;
	long size;
	float *fields;
	char* integer_or_halfinteger;
	std::ostringstream nomefile_bin, nomefile_txt;
	nomefile_bin << std::string(argv[1]);
	nomefile_txt << std::string(argv[1]) << ".txt";
	std::ifstream file_bin;
	std::ofstream file_txt;	
	file_bin.open(nomefile_bin.str().c_str(),std::ios::binary|std::ios::in);
	file_txt.open(nomefile_txt.str().c_str());
	
	if ( file_bin.fail() )
		{
			std::cout << "Input file non trovato" << std::endl;
			return -3;
		}

	file_bin.read( (char*)Ncells, 3*sizeof(int) );
	file_bin.read( (char*)rNproc, 3*sizeof(int) );
	Nproc=rNproc[0]*rNproc[1]*rNproc[2];
	xiCoords=new double[Ncells[0]];
	yiCoords=new double[Ncells[1]];
	ziCoords=new double[Ncells[2]];
	xhCoords=new double[Ncells[0]];
	yhCoords=new double[Ncells[1]];
	zhCoords=new double[Ncells[2]];
	riCoords[0]=xiCoords;
	riCoords[1]=yiCoords;
	riCoords[2]=ziCoords;
	rhCoords[0]=xhCoords;
	rhCoords[1]=yhCoords;
	rhCoords[2]=zhCoords;

	file_bin.read((char*) (&Ncomp), sizeof(int));
	
	integer_or_halfinteger = new char[3*Ncomp];	
	for (int c=0; c<Ncomp; c++){
		file_bin.read(integer_or_halfinteger+c*3,3);
	}

	for(int c=0;c<3;c++)
		file_bin.read( (char*)riCoords[c], Ncells[c]*sizeof(double));
	for(int c=0;c<3;c++)
		file_bin.read( (char*)rhCoords[c], Ncells[c]*sizeof(double));


	
	std::cout<< " Ncells:  " << Ncells[0]<<"  "<<Ncells[1]<<"  "<<Ncells[2]<<"\n";
	std::cout<< " Nprocs:  " << rNproc[0]<<"  "<<rNproc[1]<<"  "<<rNproc[2]<<"\n";
	
	/*for(int c=0; c<3;c++){
		for(int i=0;i<Ncells[c];i++)
		std::cout << riCoords[c][i]<<"  ";
		std::cout<<endl;
		}
		for(int c=0; c<3;c++){
		for(int i=0;i<Ncells[c];i++)
		std::cout << rhCoords[c][i]<<"  ";
		std::cout<<endl;
		}*/

	std::cout << "Ncomp: " << Ncomp << std::endl;
	for (int c=0; c<Ncomp; c++){
		std::cout << c << ": "
							<< (int)integer_or_halfinteger[c*3+0] << " " 
							<< (int)integer_or_halfinteger[c*3+1] << " "
							<< (int)integer_or_halfinteger[c*3+2] << std::endl;
	}

	size=((long)Ncomp)*((long)Ncells[0])*((long)Ncells[1])*((long)Ncells[2]);
	fields=new float[size];
	
	std::cout << "Reading ..." << std::endl; std::cout.flush();

	for(int rank=0;rank<Nproc;rank++){
		file_bin.read( (char*)locOrigin, 3*sizeof(int) );
		file_bin.read( (char*)locNcells, 3*sizeof(int) );
		std::cout << "rank= " << rank <<"  ";
		std::cout<<endl;
		std::cout<< "locNcells: " << locNcells[0]<<"  "<<locNcells[1]<<"  "<<locNcells[2]<<"\n";
		std::cout<< "orign: " << locOrigin[0]<<"  "<<locOrigin[1]<<"  "<<locOrigin[2]<<"\n";
		float *locFields;
		int locSize=Ncomp*locNcells[0]*locNcells[1]*locNcells[2];
		locFields=new float[locSize];
		file_bin.read( (char*)locFields, locSize*sizeof(float) );
		
		drawLoadBar(rank+1, Nproc, Nproc, 30);
	
		for(int k=0;k<locNcells[2];k++){
			for(int j=0;j<locNcells[1];j++){
				for(int i=0;i<locNcells[0];i++){
					for(int c=0;c<Ncomp;c++){
						long ii=i+locOrigin[0];
						long jj=j+locOrigin[1];
						long kk=k+locOrigin[2];
						long index    = c + Ncomp*ii + Ncomp*Ncells[0]*jj   + Ncomp*Ncells[0]*Ncells[1]*kk;
						long locIndex = c + Ncomp*i  + Ncomp*locNcells[0]*j + Ncomp*locNcells[0]*locNcells[1]*k;
						fields[index]=locFields[locIndex];
					}
				}
			}
		}
	}

	std::cout << std::endl  << "Writing to file ..." << std::endl; std::cout.flush();
	
	long k=Ncells[2]/2;

	long totPts = Ncells[0]*Ncells[1];

	for(long j=0;j<Ncells[1];j++){
		for(long i=0;i<Ncells[0];i++){

			drawLoadBar(i + j*Ncells[0] + 1, totPts, Ncells[1], 30);
			
			file_txt << setw(12) << setprecision(5) << xiCoords[i];
			file_txt << setw(12) << setprecision(5) << yiCoords[j];
			file_txt << setw(12) << setprecision(5) << ziCoords[k];
			for(int c=0;c<Ncomp;c++){
				long index = c + Ncomp*i + Ncomp*Ncells[0]*j   + Ncomp*Ncells[0]*Ncells[1]*k;
				file_txt << setw(12) << setprecision(5) << fields[index];
			}
			file_txt<<endl;
		}
	}
		
	

	file_bin.close();
	file_txt.close();

}
