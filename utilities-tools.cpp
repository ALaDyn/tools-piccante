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


#include "utilities-tools.h"

void swap_endian(float* in_f, uint64_t n)
{
  size_t i;
  union {int imio; float fmio; char arr[4];}x;
  char buff;
  for(i=0;i<n;i++)
  {
    x.fmio=in_f[i];
    buff=x.arr[0];
    x.arr[0]=x.arr[3];
    x.arr[3]=buff;
    buff=x.arr[1];
    x.arr[1]=x.arr[2];
    x.arr[2]=buff;
    in_f[i]=x.fmio;
  }
}
void swap_endian(int* in_i, int n)
{
  int i;
  union { int imio; float fmio; char arr[4]; }x;
  char buff;
  for (i = 0; i < n; i++)
  {
    x.imio = in_i[i];
    buff = x.arr[0];
    x.arr[0] = x.arr[3];
    x.arr[3] = buff;
    buff = x.arr[1];
    x.arr[1] = x.arr[2];
    x.arr[2] = buff;
    in_i[i] = x.imio;
  }
}

void swap_endian(double* in_d, int n){
  if(is_big_endian())
    return;
  int i;
  union {double frep; char arr[8];}x;
  char buff;
  for(i=0;i<n;i++){
    x.frep=in_d[i];
    buff=x.arr[0];
    x.arr[0]=x.arr[7];
    x.arr[7]=buff;

    buff=x.arr[1];
    x.arr[1]=x.arr[6];
    x.arr[6]=buff;

    buff=x.arr[2];
    x.arr[2]=x.arr[5];
    x.arr[5]=buff;

    buff=x.arr[3];
    x.arr[3]=x.arr[4];
    x.arr[4]=buff;

    in_d[i]=x.frep;
  }
}

int is_big_endian(){
  union {
    uint32_t i;
    char c[4];
  } bint = { 0x01020304 };

  return bint.c[0] == 1;
}
