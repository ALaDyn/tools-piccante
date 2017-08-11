
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

//#ifndef _XOPEN_SOURCE
//#define _XOPEN_SOURCE 500
//#endif

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <istream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <vector>
#include <random>

#if defined(CINECA)
#include <inttypes.h>
#else
#include <cstdint>
#endif

#if defined (_MSC_VER)
#include <wchar.h>
#define fseeko _fseeki64
#define ftello _ftelli64
#endif

#if defined (__MINGW32__)
#define fseeko fseeko64
#define ftello ftello64
#endif

#if defined (__CYGWIN__)
#define fseeko fseek
#define ftello ftell
#endif

using namespace std;


#define Dfloat                      float

#define MAX(x,y)                    ((x)>(y)?(x):(y))
#define MIN(x,y)                    ((x)<(y)?(x):(y))
#define THETA                       7
#define GAMMA1                      8
#define ENERGIA                     9
#define VX                          10
#define VZ                          11
#define TRUE                        1
#define FALSE                       0
//#define COMPONENTI                7               // colonne da leggere
#define NUM_COMPONENTS              7
#define NUM_QUANTITIES              10
#define VERY_BIG_POS_NUM            +1.0e30
#define VERY_BIG_NEG_NUM            -1.0e30
#define INCREASE_PLOTEXTREMS_FACTOR 0.05
#define DIMENSIONI                  12              // dimensioni lette piu' theta, gamma e energia
#define m_electron                  9.1095e-31      // electron mass
#define m_proton                    1.6726231e-27   // proton mass
//#define P_MASS                    938.272
#define m_proton_mev                938.272
#define speed_of_light              299792458
#define J2MeV                       (1e-6/1.60218e-19)


