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

#ifndef UTILITIESTOOLS_H
#define UTILITIESTOOLS_H

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<cstring>

#if defined(CINECA)
#include <inttypes.h>
#else
#include <cstdint>
#endif

void swap_endian(float* in_f, uint64_t n);
void swap_endian(int* in_i, int n);
void swap_endian(double* in_d, int n);
int is_big_endian();


#endif // UTILITIESTOOLS_H
