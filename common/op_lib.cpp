
#include <stdlib.h>                                                         
#include <stdio.h>                                                          
#include <string.h>                                                         
#include <math.h>                                                           

#ifdef _OPENMP
#include <omp.h>
#endif

#include "op_datatypes.h"

#define OP_x86

void cutilDeviceInit(int argc, char **argv){
  printf(" x86 OP2 execution \n");
}

void __syncthreads(){}

void op_fetch_data(op_dat data) {}

template <class T>
void mvHostToDevice(T **map, int size) {}

void op_decl_const_char(int, char const*, int, char*, char const*){}


#include "op_lib.cu"
