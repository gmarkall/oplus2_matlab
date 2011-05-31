
#include <stdlib.h>                                                         
#include <stdio.h>                                                          
#include <string.h>                                                         
#include <math.h>                                                           

#ifdef _OPENMP
#include <omp.h>
#endif

#include "op_lib.h"


// routines called by user code and kernels

void op_init(int argc, char **argv, int diags){
  op_init_core(argc, argv, diags);
}

void op_decl_set(int size, op_set &set, char const *name){
  op_decl_set_core(size, &set, name); 
}

void op_decl_map(op_set from, op_set to, int dim, int *map,
                 op_map &mapping, char const *name) {
  op_decl_map_core(from, to, dim, map, &mapping, name);
}

void op_fetch_data(op_dat data) {}
void op_decl_const_char(int, char const*, int, char*, char const*){}

void op_exit(){
  op_exit_core();
}

void __syncthreads(){}  // needed now in x86 kernels; get rid of it later

// routines called by op_lib_core.cpp

extern "C"
void op_mvHostToDevice(void **map, int size) {}

extern "C"
void op_cpHostToDevice(void **dat_h, void **dat_d, int size) {}


