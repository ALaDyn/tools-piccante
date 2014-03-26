#define _USE_MATH_DEFINES

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<cstring>

#if defined(CINECA)
#include <inttypes.h>
#else
#include <cstdint>
#endif

#if (defined _WIN32) || (defined _WIN64)
#include<wchar.h>
#define fseeko _fseeki64
#define ftello _ftelli64
#endif

#pragma warning(disable : 4996)


using namespace std;


#define Dfloat float 


#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define P_MASS 938.272
#define THETA 7
#define GAMMA1 8
#define ENERGIA 9
#define VX      10
#define VZ      11
#define TRUE 1
#define FALSE 0
#define COMPONENTI 7  //da leggere
#define DIMENSIONI 12 // dimensioni lette piu'theta gamma e energia


void swap_endian_f(float* , int );
int controllo(Dfloat ,Dfloat , Dfloat );
char check_domain(Dfloat , Dfloat );

