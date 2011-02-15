/*
  Open source copyright declaration based on BSD open source template:
  http://www.opensource.org/licenses/bsd-license.php

* Copyright (c) 2009, Mike Giles
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * The name of Mike Giles may not be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY Mike Giles ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL Mike Giles BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//
// header files
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "op_datatypes.h"


//
// global variables
//

int OP_set_index=0,
    OP_ptr_index=0,
    OP_dat_index=0,
    OP_nplans   =0;

op_set  * OP_set_list[10];
op_ptr  * OP_ptr_list[10];
op_dat  * OP_dat_list[10];
op_plan   OP_plans[100];


//
// run-time type-checking routine
//

int type_error(float *, op_datatype type) {
  return (type != OP_FLOAT);
}

int type_error(double *, op_datatype type) {
  return (type != OP_DOUBLE);
}

int type_error(int *, op_datatype type) {
  return (type != OP_INT);
}


//
// OP functions
//


void op_init(int argc, char **argv){
}

void op_decl_set(int size, op_set &set, char const *name){
  set.size = size;
  set.name = name;

  set.index = OP_set_index;
  OP_set_list[OP_set_index++] = &set;
}


void op_decl_ptr(op_set from, op_set to, int dim, int *ptr, op_ptr &pointer, char const *name){
  pointer.from = from;
  pointer.to   = to;
  pointer.dim  = dim;
  pointer.ptr  = ptr;
  pointer.name = name;

  pointer.index = OP_ptr_index;
  OP_ptr_list[OP_ptr_index++] = &pointer;
}


template <class T>
void op_decl_dat(op_set set, int dim, op_datatype type, T *dat, op_dat &data, char const *name){
  data.set  = set;
  data.dim  = dim;

  if (type_error(dat,type)) {
    printf("incorrect type specified for dataset \"%s\" \n",name);  exit(1);
  }

  if (type == OP_FLOAT) {
    data.fdat = (float *) dat;
  }
  else if (type == OP_DOUBLE) {
    data.ddat = (double *) dat;
  }
  else if (type == OP_INT) {
    data.idat = (int *) dat;
  }

  data.name = name;
  data.type = type;

  data.index = OP_dat_index;
  OP_dat_list[OP_dat_index++] = &data;
}

void op_decl_ddat(op_set set, int dim, op_datatype type, double *dat, op_dat &data, char const *name){
  op_decl_dat(set, dim, type, dat, data, name);
}


void op_decl_fdat(op_set set, int dim, op_datatype type, float *dat, op_dat &data, char const *name){
  op_decl_dat(set, dim, type, dat, data, name);
}


void op_decl_idat(op_set set, int dim, op_datatype type, int *dat, op_dat &data, char const *name){
  op_decl_dat(set, dim, type, dat, data, name);
}


template <class T>
void op_decl_const(int dim, op_datatype type, T *dat, char const *name){
  if (type_error(dat,type)) {
    printf("incorrect type specified for constant \"%s\" \n",name);  exit(1);
  }
}

void op_decl_dconst(int dim, op_datatype type, double *dat, char const *name){
     op_decl_const(dim, type, dat, name);
}

void op_decl_fconst(int dim, op_datatype type, float *dat, char const *name){
     op_decl_const(dim, type, dat, name);
}

void op_decl_iconst(int dim, op_datatype type, int *dat, char const *name){
     op_decl_const(dim, type, dat, name);
}


void op_fetch_data(op_dat dat){
}


//
// diagnostic output routine
//

void op_diagnostic_output(){
  if (OP_DIAGS > 1) {
    printf("\n  OP diagnostic output\n");
    printf(  "  --------------------\n");

    printf("\n       set       size\n");
    printf(  "  -------------------\n");
    for(int n=0; n<OP_set_index; n++) {
      op_set set=*OP_set_list[n];
      printf("%10s %10d\n",set.name,set.size);
    }

    printf("\n       ptr        dim       from         to\n");
    printf(  "  -----------------------------------------\n");
    for(int n=0; n<OP_ptr_index; n++) {
      op_ptr ptr=*OP_ptr_list[n];
      printf("%10s %10d %10s %10s\n",ptr.name,ptr.dim,ptr.from.name,ptr.to.name);
    }

    printf("\n       dat        dim        set\n");
    printf(  "  ------------------------------\n");
    for(int n=0; n<OP_dat_index; n++) {
      op_dat dat=*OP_dat_list[n];
      printf("%10s %10d %10s\n",dat.name,dat.dim,dat.set.name);
    }
    printf("\n");
  }
}


