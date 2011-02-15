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
#include <cutil_inline.h>
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

int   OP_consts_bytes=0;
char *OP_consts_h, *OP_consts_d;


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

extern void op_init(int argc, char **argv){
  cutilDeviceInit(argc, argv);
}

extern void op_decl_set(int size, op_set &set, char const *name){
  set.size = size;
  set.name = name;

  set.index = OP_set_index;
  OP_set_list[OP_set_index++] = &set;
}

extern void op_decl_ptr(op_set from, op_set to, int dim, int *ptr, op_ptr &pointer, char const *name){
  pointer.from = from;
  pointer.to   = to;
  pointer.dim  = dim;
  pointer.ptr  = ptr;
  pointer.name = name;

  pointer.index = OP_ptr_index;
  OP_ptr_list[OP_ptr_index++] = &pointer;
}


template <class T>
void op_decl_dat_T(op_set set, int dim, op_datatype type, T *dat, op_dat &data, char const *name){

  if (type_error(dat,type)) {
    printf("incorrect type specified for dataset \"%s\" \n",name);  exit(1);
  }

  data.set   = set;
  data.dim   = dim;
  data.dat   = (char *) dat;
  data.name  = name;
  data.type  = type;
  data.size  = dim*sizeof(T);
  data.index = OP_dat_index;
  OP_dat_list[OP_dat_index++] = &data;

  cutilSafeCall(cudaMalloc((void **)&data.dat_d, data.size*set.size));
  cutilSafeCall(cudaMemcpy(data.dat_d, data.dat, data.size*set.size,
                cudaMemcpyHostToDevice));
}

extern
void op_decl_dat(op_set set, int dim, op_datatype type, double *dat, op_dat &data, char const *name){
  op_decl_dat_T(set, dim, type, dat, data, name);
}

extern
void op_decl_dat(op_set set, int dim, op_datatype type, float *dat, op_dat &data, char const *name){
  op_decl_dat_T(set, dim, type, dat, data, name);
}

extern
void op_decl_dat(op_set set, int dim, op_datatype type, int *dat, op_dat &data, char const *name){
  op_decl_dat_T(set, dim, type, dat, data, name);
}


template <class T>
void op_decl_const_T(int dim, op_datatype type, T *dat, char const *name){
  if (type_error(dat,type)) {
    printf("incorrect type specified for constant \"%s\" \n",name);  exit(1);
  }

  // printf(" op_decl_const: name = %s, size = %d\n",name,sizeof(T)*dim);
  cutilSafeCall( cudaMemcpyToSymbol(name, dat, sizeof(T)*dim));
}

extern void op_decl_const(int dim, op_datatype type, double *dat, char const *name){
            op_decl_const_T(dim, type, dat, name);
}

extern void op_decl_const(int dim, op_datatype type, float *dat, char const *name){
            op_decl_const_T(dim, type, dat, name);
}

extern void op_decl_const(int dim, op_datatype type, int *dat, char const *name){
            op_decl_const_T(dim, type, dat, name);
}

//
// diagnostic output routine
//

extern void op_diagnostic_output(){
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


// comparison function for integer quicksort

int comp(const void *a2, const void *b2) {
  int *a = (int *)a2;
  int *b = (int *)b2;

  if (*a == *b)
    return 0;
  else
    if (*a < *b)
      return -1;
    else
      return 1;
}

//
// utility routine to move arrays to GPU device
//

template <class T>
void mvHostToDevice(T **ptr, int size) {
  T *tmp;
  cutilSafeCall(cudaMalloc((void **)&tmp, size));
  cutilSafeCall(cudaMemcpy(tmp, *ptr, size, cudaMemcpyHostToDevice));
  free(*ptr);
  *ptr = tmp;
}


//
// utility routine to copy dataset back to host
//

extern void op_fetch_data(op_dat data) {
  cutilSafeCall(cudaMemcpy(data.dat, data.dat_d, data.size*data.set.size,
                cudaMemcpyDeviceToHost));
  cutilSafeCall(cudaThreadSynchronize());
}


//
// utility routine to resize constant arrays, if necessary
//

void reallocConstArrays(int consts_bytes) {
  if (OP_consts_bytes>0) {
    free(OP_consts_h);
    cutilSafeCall(cudaFree(OP_consts_d));
  }
  OP_consts_bytes = 4*consts_bytes;
  OP_consts_h = (char *) malloc(OP_consts_bytes);
  cutilSafeCall(cudaMalloc((void **)&OP_consts_d, OP_consts_bytes));
}

//
// utility routine to move constant arrays, if necessary
//

void mvConstArraysToDevice(int consts_bytes) {
  cutilSafeCall(cudaMemcpy(OP_consts_d, OP_consts_h, consts_bytes,
                cudaMemcpyHostToDevice));
}


//
// find existing execution plan, or construct a new one
//

extern op_plan * plan(char const * name, op_set set, int nargs, op_dat *args, int *idxs,
      op_ptr *ptrs, int *dims, op_datatype *typs, op_access *accs, int ninds, int *inds){

  // first look for an existing execution plan

  int ip=0, match=0;

  while (match==0 && ip<OP_nplans) {
    if ( (strcmp(name,        OP_plans[ip].name)==0)
             && (set.index == OP_plans[ip].set_index)
             && (nargs     == OP_plans[ip].nargs) ) {
      match = 1;
      for (int m=0; m<nargs; m++) {
        match = match && (args[m].index == OP_plans[ip].arg_idxs[m])
                      && (idxs[m]       == OP_plans[ip].idxs[m])
                      && (ptrs[m].index == OP_plans[ip].ptr_idxs[m])
                      && (dims[m]       == OP_plans[ip].dims[m])
                      && (typs[m]       == OP_plans[ip].typs[m])
                      && (accs[m]       == OP_plans[ip].accs[m]);
      }
    }
    ip++;
  }

  if (match) {
    ip--;
    if (OP_DIAGS > 1) printf(" old execution plan #%d\n",ip);
    return &(OP_plans[ip]);
  }
  else {
    if (OP_DIAGS > 1) printf(" new execution plan #%d\n",ip);
  }

  // consistency checks

  if (OP_DIAGS > 0) {
    for (int m=0; m<nargs; m++) {
      if (idxs[m] == -1) {
        //if (ptrs[m].index != -1) {
        if (ptrs[m].ptr != NULL) {
          printf("error2: wrong pointer for arg %d in kernel \"%s\"\n",m,name);
          printf("ptrs[m].index = %d\n",ptrs[m].index);
          printf("ptrs[m].name  = %s\n",ptrs[m].name);
          exit(1);
        }
      }
      else {
        if (set.index         != ptrs[m].from.index ||
            args[m].set.index != ptrs[m].to.index) {
          printf("error: wrong pointer for arg %d in kernel \"%s\"\n",m,name);
          exit(1);
        }
        if (ptrs[m].dim <= idxs[m]) {
          printf(" %d %d",ptrs[m].dim,idxs[m]);
          printf("error: invalid pointer index for arg %d in kernel \"%s\"\n",m,name);
          exit(1);
        }
      }
      if (args[m].type != typs[m]) {
        printf("error: wrong datatype for arg %d in kernel \"%s\"\n",m,name);
        exit(1);
      }
      if (args[m].dim != dims[m] && args[m].set.size>0) {
        printf("error: wrong dimension for arg %d in kernel \"%s\"\n",m,name);
        exit(1);
      }
    }
  }

  // set blocksize and number of blocks

  // int bsize   = 1024;   // blocksize
  int bsize   = 100;   // blocksize
  int nblocks = (set.size-1)/bsize + 1;

  printf(" number of blocks = %d\n",nblocks);

  // allocate memory for new execution plan and store input arguments

  OP_plans[ip].arg_idxs  = (int *)malloc(nargs*sizeof(int));
  OP_plans[ip].idxs      = (int *)malloc(nargs*sizeof(int));
  OP_plans[ip].ptr_idxs  = (int *)malloc(nargs*sizeof(int));
  OP_plans[ip].dims      = (int *)malloc(nargs*sizeof(int));
  OP_plans[ip].typs      = (op_datatype *)malloc(nargs*sizeof(op_datatype));
  OP_plans[ip].accs      = (op_access *)malloc(nargs*sizeof(op_access));

  OP_plans[ip].nthrcol   = (int *)malloc(nblocks*sizeof(int));
  OP_plans[ip].thrcol    = (int *)calloc(set.size,sizeof(int));
  OP_plans[ip].offset    = (int *)malloc(nblocks*sizeof(int));
  OP_plans[ip].ind_ptrs  = (int **)malloc(ninds*sizeof(int *));
  OP_plans[ip].ind_offs  = (int **)malloc(ninds*sizeof(int *));
  OP_plans[ip].ind_sizes = (int **)malloc(ninds*sizeof(int *));
  OP_plans[ip].ptrs      = (int **)malloc(nargs*sizeof(int *));
  OP_plans[ip].nelems    = (int *)malloc(nblocks*sizeof(int));
  OP_plans[ip].ncolblk   = (int *)calloc(set.size,sizeof(int)); // max possibly needed
  OP_plans[ip].blkmap    = (int *)calloc(nblocks,sizeof(int));

  for (int m=0; m<ninds; m++) {
    int count = 0;
    for (int m2=0; m2<nargs; m2++)
      if (inds[m2]==m) count++;
    OP_plans[ip].ind_ptrs[m]  = (int *)malloc(count*set.size*sizeof(int));
    OP_plans[ip].ind_offs[m]  = (int *)malloc(nblocks*sizeof(int));
    OP_plans[ip].ind_sizes[m] = (int *)malloc(nblocks*sizeof(int));
  }

  for (int m=0; m<nargs; m++) {
    OP_plans[ip].ptrs[m]     = (int *)malloc(set.size*sizeof(int));

    OP_plans[ip].arg_idxs[m] = args[m].index;
    OP_plans[ip].idxs[m]     = idxs[m];
    OP_plans[ip].ptr_idxs[m] = ptrs[m].index;
    OP_plans[ip].dims[m]     = dims[m];
    OP_plans[ip].typs[m]     = typs[m];
    OP_plans[ip].accs[m]     = accs[m];
  }

  OP_plans[ip].name      = name;
  OP_plans[ip].set_index = set.index;
  OP_plans[ip].nargs     = nargs;
    
  OP_nplans++;


  // allocate working arrays

  uint **work;
  work = (uint **)malloc(ninds*sizeof(uint *));

  for (int m=0; m<ninds; m++) {
    int m2 = 0;
    while(inds[m2]!=m) m2++;

    work[m] = (uint *)malloc(ptrs[m2].to.size*sizeof(uint));
  }

  //  uint **work;
  //  work = (uint **)malloc(nargs*sizeof(uint *));
  //
  //  for (int m=0; m<nargs; m++) {
  //    work[m] = (uint *)malloc(ptrs[m].to.size*sizeof(uint));
  //  }

  int *work2;
  work2 = (int *)malloc(nargs*bsize*sizeof(int));  // max possibly needed


  // process set one block at a time

  int *nindirect;
  nindirect = (int *)calloc(ninds,sizeof(int));  // total number of indirect elements

  for (int b=0; b<nblocks; b++) {
    int bs = MIN(bsize, set.size - b*bsize);

    OP_plans[ip].offset[b] = b*bsize;    // offset for block
    OP_plans[ip].nelems[b] = bs;         // size of block

    // loop over indirection sets

    for (int m=0; m<ninds; m++) {

      // build the list of elements indirectly referenced in this block

      int ne = 0;  // number of elements
      for (int m2=0; m2<nargs; m2++) {
        if (inds[m2]==m) {
          for (int e=b*bsize; e<b*bsize+bs; e++)
            work2[ne++] = ptrs[m2].ptr[e];
	}
      }

      // sort them, then eliminate duplicates

      qsort(work2,ne,sizeof(int),comp);
        
      int e = 0;
      int p = 0;
      while (p<ne) {
        work2[e] = work2[p];
        while (p<ne && work2[p]==work2[e]) p++;
        e++;
      }
      ne = e;  // number of distinct elements

      if (OP_DIAGS > 5) {
        printf(" indirection set %d: ",m);
        for (int e=0; e<ne; e++) printf(" %d",work2[e]);
        printf(" \n");
      }

      // store mapping and renumbered pointers in execution plan

      for (int e=0; e<ne; e++) {
        OP_plans[ip].ind_ptrs[m][nindirect[m]++] = work2[e];
        work[m][work2[e]] = e;   // inverse mapping
      }

      for (int m2=0; m2<nargs; m2++) {
        if (inds[m2]==m) {
          for (int e=b*bsize; e<b*bsize+bs; e++)
            OP_plans[ip].ptrs[m2][e] = work[m][ptrs[m2].ptr[e]];
	}
      }

      if (b==0) {
        OP_plans[ip].ind_offs[m][b]  = 0;
        OP_plans[ip].ind_sizes[m][b] = nindirect[m];
      }
      else {
        OP_plans[ip].ind_offs[m][b]  = OP_plans[ip].ind_offs[m][b-1]
                                     + OP_plans[ip].ind_sizes[m][b-1];
        OP_plans[ip].ind_sizes[m][b] = nindirect[m] - OP_plans[ip].ind_offs[m][b];
      }

      // printf("b, m, nindirect[m] = %d %d %d\n",b,m,nindirect[m]);
    }


    // print out re-numbered pointers

    /*
    for (int m=0; m<nargs; m++) {
      if (inds[m]>=0) {
        printf(" pointer table %d\n",m);
        for (int e=0; e<set.size; e++)
          printf(" ptr = %d\n",OP_plans[ip].ptrs[m][e]);
      }
    }

    for (int m=0; m<ninds; m++) {
      printf(" indirect set %d\n",m);
      for (int b=0; b<nblocks; b++) {
        printf("OP_plans[ip].ind_sizes[m][b] = %d\n", OP_plans[ip].ind_sizes[m][b]);
        printf("OP_plans[ip].ind_offs[m][b] = %d\n", OP_plans[ip].ind_offs[m][b]);
      }
    }
    */

    // now colour main set elements

    int repeat  = 1;
    int ncolor  = 0;
    int ncolors = 0;

    while (repeat) {
      repeat = 0;

      for (int m=0; m<nargs; m++) {
        if (inds[m]>=0)
          for (int e=b*bsize; e<b*bsize+bs; e++)
            work[inds[m]][ptrs[m].ptr[e]] = 0;  // zero out relevant bits of color arrays
      }

      for (int e=b*bsize; e<b*bsize+bs; e++) {
        if (OP_plans[ip].thrcol[e]==0) {
          int mask = 0;
          for (int m=0; m<nargs; m++)
            if (inds[m]>=0 && accs[m]==OP_INC)
              mask |= work[inds[m]][ptrs[m].ptr[e]];    // set bits of mask

          int color = ffs(~mask) - 1;   // find first bit not set
          // printf("block, thread, color = %d, %d, %d \n",b,e,color);
          if (color==-1) {              // run out of colors on this pass
            repeat = 1;
          }
          else {
            OP_plans[ip].thrcol[e] = ncolor+color;
            mask    = 1 << color;
            ncolors = MAX(ncolors, ncolor+color+1);

            for (int m=0; m<nargs; m++)
              if (inds[m]>=0 && accs[m]==OP_INC)
                work[inds[m]][ptrs[m].ptr[e]] |= mask;  // set bit of color array
          }
        }
      }

      ncolor += 32;   // increment base level
    }

    OP_plans[ip].nthrcol[b] = ncolors;  // number of thread colors in this block

    // reorder elements by color?

  }


  // colour the blocks, after initialising colors to 0

  int *blk_col;
  blk_col = (int *)calloc(nblocks,sizeof(int));

  int repeat  = 1;
  int ncolor  = 0;
  int ncolors = 0;

  while (repeat) {
    repeat = 0;

    for (int m=0; m<nargs; m++) {
      if (inds[m]>=0) 
        for (int e=0; e<ptrs[m].to.size; e++)
          work[inds[m]][e] = 0;               // zero out color arrays
    }

    for (int b=0; b<nblocks; b++) {
      if (blk_col[b] == 0) {          // color not yet assigned to block
        int  bs   = MIN(bsize, set.size - b*bsize);
        uint mask = 0;

        for (int m=0; m<nargs; m++) {
          if (inds[m]>=0) 
            for (int e=b*bsize; e<b*bsize+bs; e++)
              mask |= work[inds[m]][ptrs[m].ptr[e]];    // set bits of mask
        }

        int color = ffs(~mask) - 1;   // find first bit not set
        // printf("block, color = %d, %d \n",b,color);
        if (color==-1) {              // run out of colors on this pass
          repeat = 1;
        }
        else {
          blk_col[b] = ncolor + color;
          mask    = 1 << color;
          ncolors = MAX(ncolors, ncolor+color+1);

          for (int m=0; m<nargs; m++) {
            if (inds[m]>=0)
              for (int e=b*bsize; e<b*bsize+bs; e++)
                work[inds[m]][ptrs[m].ptr[e]] |= mask;
          }
        }
      }
    }

    ncolor += 32;   // increment base level
  }


  // store block mapping and number of blocks per color


  OP_plans[ip].ncolors = ncolors;

  for (int b=0; b<nblocks; b++) {
    // printf("b, blk_col = %d, %d\n",b,blk_col[b]);
    OP_plans[ip].ncolblk[blk_col[b]]++;  // number of blocks of each color
  }

  for (int c=1; c<ncolors; c++)
    OP_plans[ip].ncolblk[c] += OP_plans[ip].ncolblk[c-1]; // cumsum

  for (int c=0; c<ncolors; c++) work2[c]=0;

  for (int b=0; b<nblocks; b++) {
    int c  = blk_col[b];
    int b2 = work2[c];     // number of preceding blocks of this color
    if (c>0) b2 += OP_plans[ip].ncolblk[c-1];  // plus previous colors

    OP_plans[ip].blkmap[b2] = b;

    work2[c]++;            // increment counter
  }

  for (int c=ncolors-1; c>0; c--)
    OP_plans[ip].ncolblk[c] -= OP_plans[ip].ncolblk[c-1]; // undo cumsum

  /*
  printf("ncolblk = ");
  for (int c=0; c<ncolors; c++)
    printf("%d ",OP_plans[ip].ncolblk[c]);
  printf("\n");
  */

  // reorder blocks by color?


  // work out shared memory requirements

  // printf(" working out shared memory needs\n");

  OP_plans[ip].nshared = 0;

  for (int b=0; b<nblocks; b++) {
    int nbytes = 0;
    for (int m=0; m<ninds; m++) {
      int m2 = 0;
      while(inds[m2]!=m) m2++;

      int nelems = OP_plans[ip].ind_sizes[m][b]*dims[m2];
      if (typs[m2]==OP_FLOAT)
        nbytes += nelems*sizeof(float);
      else if (typs[m2]==OP_DOUBLE)
        nbytes += nelems*sizeof(double);
      else if (typs[m2]==OP_INT)
        nbytes += nelems*sizeof(int);
    }
    OP_plans[ip].nshared = MAX(OP_plans[ip].nshared,nbytes);
  }

  // printf(" shared memory = %d bytes \n",OP_plans[ip].nshared);

  // move plan arrays to GPU

  for (int m=0; m<ninds; m++) {
    mvHostToDevice(&(OP_plans[ip].ind_ptrs[m]), sizeof(int)*nindirect[m]);
    mvHostToDevice(&(OP_plans[ip].ind_sizes[m]),sizeof(int)*nblocks);
    mvHostToDevice(&(OP_plans[ip].ind_offs[m]), sizeof(int)*nblocks);
  }

  for (int m=0; m<nargs; m++) {
    if (inds[m]>=0)
      mvHostToDevice(&(OP_plans[ip].ptrs[m]), sizeof(int)*set.size);
  }

  mvHostToDevice(&(OP_plans[ip].nthrcol),sizeof(int)*nblocks);
  mvHostToDevice(&(OP_plans[ip].thrcol ),sizeof(int)*set.size);
  mvHostToDevice(&(OP_plans[ip].offset ),sizeof(int)*nblocks);
  mvHostToDevice(&(OP_plans[ip].nelems ),sizeof(int)*nblocks);
  mvHostToDevice(&(OP_plans[ip].blkmap ),sizeof(int)*nblocks);

  // free work arrays

  for (int m=0; m<ninds; m++) free(work[m]);
  free(work);
  free(work2);
  free(blk_col);
  free(nindirect);

  // return pointer to plan

  return &(OP_plans[ip]);
}