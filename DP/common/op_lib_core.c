/*
  Open source copyright declaration based on BSD open source template:
  http://www.opensource.org/licenses/bsd-license.php

* Copyright (c) 2009-2011, Mike Giles
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
// OP2 core header file -- defines datatypes and core functions
//

#include "op_lib_core.h"

//
// global variables
//

int OP_diags          =0,
    OP_part_size      =0,
    OP_block_size     =512,
    OP_cache_line_size=128;

int OP_set_index =0, OP_set_max =0,
    OP_map_index =0, OP_map_max =0,
    OP_dat_index =0, OP_dat_max =0,
    OP_plan_index=0, OP_plan_max=0,
                     OP_kern_max=0;

op_set    *OP_set_list;
op_map    *OP_map_list;
op_dat    *OP_dat_list;
op_plan   *OP_plans;
op_kernel *OP_kernels;

//
// function prototypes
//

extern void op_mvHostToDevice(void **, int);
extern void op_cpHostToDevice(void **, void **, int);

//
// OP core functions
//

void op_init_core(int argc, char **argv, int diags){
  OP_diags = diags;

#ifdef OP_BLOCK_SIZE
  OP_block_size = OP_BLOCK_SIZE;
#endif
#ifdef OP_PART_SIZE
  OP_part_size = OP_PART_SIZE;
#endif

  for (int n=1; n<argc; n++) {
    if (strncmp(argv[n],"OP_BLOCK_SIZE=",14)==0) {
      OP_block_size = atoi(argv[n]+14);
      printf("\n OP_block_size = %d \n", OP_block_size);
    }
    if (strncmp(argv[n],"OP_PART_SIZE=",13)==0) {
      OP_part_size = atoi(argv[n]+13);
      printf("\n OP_part_size  = %d \n", OP_part_size);
    }
    if (strncmp(argv[n],"OP_CACHE_LINE_SIZE=",19)==0) {
      OP_cache_line_size = atoi(argv[n]+19);
      printf("\n OP_cache_line_size  = %d \n", OP_cache_line_size);
    }
  }
}

op_set op_decl_set(int size, char const *name){

  if (size<=0) {
    printf(" op_decl_set error -- negative/zero size for set: %s\n",name);
    exit(-1);
  }

  if (OP_set_index==OP_set_max) {
    OP_set_max += 10;
    OP_set_list = (op_set *)realloc(OP_set_list,
                                    OP_set_max*sizeof(op_set));
    if (OP_set_list==NULL) {
      printf(" op_decl_set error -- error reallocating memory\n");
      exit(-1);  
    }
  }

  op_set set = (op_set) malloc(sizeof(op_set_core));
  set->index = OP_set_index;
  set->size  = size;
  set->name  = name;

  OP_set_list[OP_set_index++] = set;

  return set;
}

op_map op_decl_map(op_set from, op_set to, int dim, int *imap,
                                            char const *name){
  if (from==NULL) {
    printf(" op_decl_map error -- invalid 'from' set for map %s\n",name);
    exit(-1);
  }

  if (to==NULL) {
    printf("op_decl_map error -- invalid 'to' set for map %s\n",name);
    exit(-1);
  }

  if (dim<=0) {
    printf("op_decl_map error -- negative/zero dimension for map %s\n",name);
    exit(-1);
  }

  /*for (int d=0; d<dim; d++) {
    for (int n=0; n<from->size; n++) {
      if (imap[d+n*dim]<0 || imap[d+n*dim]>=to->size) {
        printf("op_decl_map error -- invalid data for map %s\n",name);
        printf("element = %d, dimension = %d, map = %d\n",n,d,imap[d+n*dim]);
        exit(-1);
      }
    }
  }*/

  if (OP_map_index==OP_map_max) {
    OP_map_max += 10;
    OP_map_list = (op_map *) realloc(OP_map_list,
                                     OP_map_max*sizeof(op_map));
    if (OP_map_list==NULL) {
      printf(" op_decl_map error -- error reallocating memory\n");
      exit(-1);  
    }
  }

  op_map map = (op_map) malloc(sizeof(op_map_core));
  map->index = OP_map_index;
  map->from  = from;
  map->to    = to;
  map->dim   = dim;
  map->map   = imap;
  map->name  = name;

  OP_map_list[OP_map_index++] = map;

  return map;
}

op_dat op_decl_dat_core(op_set set, int dim, char const *type,
                       int size, char *data, char const *name){

  if (set==NULL) {
    printf("op_decl_dat error -- invalid set for data: %s\n",name);
    exit(-1);
  }

  if (dim<=0) {
    printf("op_decl_dat error -- negative/zero dimension for data: %s\n",name);
    exit(-1);
  }

  if (OP_dat_index==OP_dat_max) {
    OP_dat_max += 10;
    OP_dat_list = (op_dat *) realloc(OP_dat_list,
                                   OP_dat_max*sizeof(op_dat));
    if (OP_dat_list==NULL) {
      printf(" op_decl_dat error -- error reallocating memory\n");
      exit(-1);  
    }
  }

  op_dat dat  = (op_dat) malloc(sizeof(op_dat_core));
  dat->index  = OP_dat_index;
  dat->set    = set;
  dat->dim    = dim;
  dat->data   = data;
  dat->data_d = NULL;
  dat->name   = name;
  dat->type   = type;
  dat->size   = dim*size;

  OP_dat_list[OP_dat_index++] = dat;

  return dat;
}

void op_exit_core(){

  // free storage for plans

  for(int ip=0; ip<OP_plan_index; ip++) {
    free(OP_plans[ip].dats);
    free(OP_plans[ip].idxs);
    free(OP_plans[ip].maps);
    free(OP_plans[ip].accs);
    free(OP_plans[ip].ind_maps);
    free(OP_plans[ip].nindirect);
    free(OP_plans[ip].loc_maps);
    free(OP_plans[ip].ncolblk);
  }

  free(OP_plans);

  // free storage and pointers for sets, maps and data

  for(int i=0; i<OP_set_index; i++) {
    free(OP_set_list[i]);
  }
  free(OP_set_list);

  for(int i=0; i<OP_map_index; i++) {
    free(OP_map_list[i]);
  }
  free(OP_map_list);

  for(int i=0; i<OP_dat_index; i++) {
    free(OP_dat_list[i]);
  }
  free(OP_dat_list);

  // free storage for timing info

  free(OP_kernels);

  // reset initial values

  OP_set_index =0; OP_set_max =0;
  OP_map_index =0; OP_map_max =0;
  OP_dat_index =0; OP_dat_max =0;
  OP_plan_index=0; OP_plan_max=0;
  OP_kern_max=0;
}

//
// op_arg routines
//

void op_err_print(const char *error_string, int m, const char *name) {
  printf("error: arg %d in kernel \"%s\"\n",m,name);
  printf("%s \n", error_string);
  exit(1);
}

void op_arg_check(op_set set,int m,op_arg arg,int *ninds,const char *name){

  // error checking for op_arg_dat

  if (arg.argtype==OP_ARG_DAT) {
    if (set==NULL)
      op_err_print("invalid set",m,name);

    if (arg.map==NULL && arg.dat->set!=set)
      op_err_print("dataset set does not match loop set",m,name);

    if (arg.map!=NULL && (arg.map->from!=set || arg.map->to!=arg.dat->set))
      op_err_print("mapping error",m,name);

    if ( (arg.map==NULL &&  arg.idx!=-1) ||
         (arg.map!=NULL && (arg.idx<0 || arg.idx>=arg.map->dim) ) )
      op_err_print("invalid index",m,name);

    if (arg.dat->dim != arg.dim)
      op_err_print("dataset dim does not match declared dim",m,name);

    if (strcmp(arg.dat->type,arg.type))
      op_err_print("dataset type does not match declared type",m,name);

    if (arg.idx>=0) (*ninds)++;
  }

  // error checking for op_arg_gbl

  if (arg.argtype==OP_ARG_GBL) {
    if (!strcmp(arg.type,"error"))
      op_err_print("datatype does not match declared type",m,name);

    if (arg.dim<=0)
      op_err_print("dimension should be strictly positive",m,name);

    if (arg.data==NULL)
      op_err_print("NULL pointer for global data",m,name);
  }
}

op_arg op_arg_dat(op_dat dat, int idx, op_map map,
                  int dim, const char *typ, op_access acc){
  op_arg arg;

  arg.argtype = OP_ARG_DAT;

  arg.dat  = dat;
  arg.map  = map;
  arg.dim  = dim;
  arg.idx  = idx;
  if (dat!=NULL) {
    arg.size   = dat->size;
    arg.data   = dat->data;
    arg.data_d = dat->data_d;
  }
  arg.type = typ;
  arg.acc  = acc;

  return arg;
}


op_arg op_arg_gbl_core(char *data, int dim,
                       const char *typ, op_access acc){
  op_arg arg;

  arg.argtype = OP_ARG_GBL;

  arg.dat  = NULL;
  arg.map  = NULL;
  arg.dim  = dim;
  arg.idx  = -1;
  arg.size = 0;
  arg.data = data;
  arg.type = typ;
  arg.acc  = acc;

  return arg;
}

//
// diagnostic routines
//

void op_diagnostic_output(){
  if (OP_diags > 1) {
    printf("\n  OP diagnostic output\n");
    printf(  "  --------------------\n");

    printf("\n       set       size\n");
    printf(  "  -------------------\n");
    for(int n=0; n<OP_set_index; n++) {
      printf("%10s %10d\n",
             OP_set_list[n]->name,OP_set_list[n]->size);
    }

    printf("\n       map        dim       from         to\n");
    printf(  "  -----------------------------------------\n");
    for(int n=0; n<OP_map_index; n++) {
      printf("%10s %10d %10s %10s\n",
             OP_map_list[n]->name,OP_map_list[n]->dim,
             OP_map_list[n]->from->name,OP_map_list[n]->to->name);
    }

    printf("\n       dat        dim        set\n");
    printf(  "  ------------------------------\n");
    for(int n=0; n<OP_dat_index; n++) {
      printf("%10s %10d %10s\n", OP_dat_list[n]->name,
             OP_dat_list[n]->dim,OP_dat_list[n]->set->name);
    }
    printf("\n");
  }
}

void op_timing_output() {
  
  double max_time = 0.0;
  if (OP_kern_max>0) {
    printf("\n  count     time     GB/s     GB/s   kernel name ");
    printf("\n ----------------------------------------------- \n");
    for (int n=0; n<OP_kern_max; n++) {
      if (OP_kernels[n].count>0) {
        if (OP_kernels[n].transfer2==0.0f)
          printf(" %6d  %8.4f %8.4f            %s \n",
	       OP_kernels[n].count,
               OP_kernels[n].time,
               OP_kernels[n].transfer/(1e9f*OP_kernels[n].time),
               OP_kernels[n].name);
        else
          printf(" %6d  %8.4f %8.4f %8.4f   %s \n",
	       OP_kernels[n].count,
               OP_kernels[n].time,
               OP_kernels[n].transfer/(1e9f*OP_kernels[n].time),
               OP_kernels[n].transfer2/(1e9f*OP_kernels[n].time),
               OP_kernels[n].name);
          max_time = max_time + OP_kernels[n].time;
      }
    }
    printf("Max total runtime = %f\n",max_time);
  }
}

void op_timing_realloc(int kernel){
  int OP_kern_max_new; 

  if (kernel>=OP_kern_max) {
    // printf("allocating more memory for OP_kernels \n");
    OP_kern_max_new = kernel + 10;
    OP_kernels = (op_kernel *) realloc(OP_kernels,
                OP_kern_max_new*sizeof(op_kernel));
    if (OP_kernels==NULL) {
      printf(" op_timing_realloc error \n");
      exit(-1);  
    }

    for (int n=OP_kern_max; n<OP_kern_max_new; n++) {
      OP_kernels[n].count     = 0;
      OP_kernels[n].time      = 0.0f;
      OP_kernels[n].transfer  = 0.0f;
      OP_kernels[n].transfer2 = 0.0f;
      OP_kernels[n].name      = "unused";
    }
    OP_kern_max = OP_kern_max_new;
  }
}


//
// routines below are needed for generated CUDA and OpenMP code
//

//
// timing routine from Gihan Mudalige
//

#include <sys/time.h>

void op_timers(double *cpu, double *et) {
  struct timeval t;

  gettimeofday( &t, (struct timezone *)0 );
  *et = t.tv_sec + t.tv_usec*1.0e-6;
}

//
// comparison function for integer quicksort in op_plan
//

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
// plan check routine
//

void op_plan_check(op_plan OP_plan, int ninds, int *inds) {

  int err, ntot;

  int nblock = 0;
  for (int col=0; col<OP_plan.ncolors; col++) nblock += OP_plan.ncolblk[col];

  //
  // check total size
  //

  int nelem = 0;
  for (int n=0; n<nblock; n++) nelem += OP_plan.nelems[n];

  if (nelem != OP_plan.set->size) {
    printf(" *** OP_plan_check: nelems error \n");
  }
  else if (OP_diags>6) {
    printf(" *** OP_plan_check: nelems   OK \n");
  }

  //
  // check offset and nelems are consistent
  //

  err  = 0;
  ntot = 0;

  for (int n=0; n<nblock; n++) {
    err  += (OP_plan.offset[n] != ntot);
    ntot +=  OP_plan.nelems[n];
  }

  if (err != 0) {
    printf(" *** OP_plan_check: offset error \n");
  }
  else if (OP_diags>6) {
    printf(" *** OP_plan_check: offset   OK \n");
  }

  //
  // check blkmap permutation
  //

  int *blkmap = (int *) malloc(nblock*sizeof(int));
  for (int n=0; n<nblock; n++) blkmap[n] = OP_plan.blkmap[n];
  qsort(blkmap,nblock,sizeof(int),comp);

  err = 0;
  for (int n=0; n<nblock; n++) err += (blkmap[n] != n);

  free(blkmap);

  if (err != 0) {
    printf(" *** OP_plan_check: blkmap error \n");
  }
  else if (OP_diags>6) {
    printf(" *** OP_plan_check: blkmap   OK \n");
  }

  //
  // check ind_offs and ind_sizes are consistent
  //

  err  = 0;

  for (int i = 0; i<ninds; i++) {
    ntot = 0;

    for (int n=0; n<nblock; n++) {
      err  += (OP_plan.ind_offs[i+n*ninds] != ntot);
      ntot +=  OP_plan.ind_sizes[i+n*ninds];
    }
  }

  if (err != 0) {
    printf(" *** OP_plan_check: ind_offs error \n");
  }
  else if (OP_diags>6) {
    printf(" *** OP_plan_check: ind_offs OK \n");
  }

  //
  // check ind_maps correctly ordered within each block
  // and indices within range
  //

  err = 0;

  for (int m = 0; m<ninds; m++) {
    int m2 = 0;
    while(inds[m2]!=m) m2++;
    int set_size = OP_plan.maps[m2]->to->size;

    ntot = 0;

    for (int n=0; n<nblock; n++) {
      int last = -1;
      for (int e=ntot; e<ntot+OP_plan.ind_sizes[m+n*ninds]; e++) {
        err  += (OP_plan.ind_maps[m][e] <= last);
        last  = OP_plan.ind_maps[m][e]; 
      }
      err  += (last >= set_size);
      ntot +=  OP_plan.ind_sizes[m+n*ninds];
    }
  }

  if (err != 0) {
    printf(" *** OP_plan_check: ind_maps error \n");
  }
  else if (OP_diags>6) {
    printf(" *** OP_plan_check: ind_maps OK \n");
  }

  //
  // check maps (most likely source of errors)
  //

  err = 0;

  for (int m=0; m<OP_plan.nargs; m++) {
    if (OP_plan.maps[m]!=NULL) {
      op_map map = OP_plan.maps[m];
      int    m2  = inds[m];

      ntot = 0;
      for (int n=0; n<nblock; n++) {
        for (int e=ntot; e<ntot+OP_plan.nelems[n]; e++) {
          int p_local  = OP_plan.loc_maps[m][e];
          int p_global = OP_plan.ind_maps[m2][p_local+OP_plan.ind_offs[m2+n*ninds]]; 

          err += (p_global != map->map[OP_plan.idxs[m]+e*map->dim]);
        }
        ntot +=  OP_plan.nelems[n];
      }
    }
  }

  if (err != 0) {
    printf(" *** OP_plan_check: maps error \n");
  }
  else if (OP_diags>6) {
    printf(" *** OP_plan_check: maps     OK \n");
  }


  //
  // check thread and block coloring
  //

  return;
}


//
// OP plan construction
//

op_plan * op_plan_core(char const *name, op_set set, int part_size,
                     int nargs, op_arg *args, int ninds, int *inds){

  // first look for an existing execution plan

  int ip=0, match=0;

  while (match==0 && ip<OP_plan_index) {
    if ( (strcmp(name,        OP_plans[ip].name) == 0)
             && (set       == OP_plans[ip].set       )
             && (nargs     == OP_plans[ip].nargs     )
             && (ninds     == OP_plans[ip].ninds     )
             && (part_size == OP_plans[ip].part_size ) ) {
      match = 1;
      for (int m=0; m<nargs; m++) {
        match = match && (args[m].dat == OP_plans[ip].dats[m])
                      && (args[m].map == OP_plans[ip].maps[m])
                      && (args[m].idx == OP_plans[ip].idxs[m])
                      && (args[m].acc == OP_plans[ip].accs[m]);
      }
    }
    ip++;
  }

  if (match) {
    ip--;
    if (OP_diags > 3) printf(" old execution plan #%d\n",ip);
    OP_plans[ip].count++;
    return &(OP_plans[ip]);
  }
  else {
    if (OP_diags > 1) printf(" new execution plan #%d for kernel %s\n",ip,name);
  }

  // work out worst case shared memory requirement per element

  int maxbytes = 0;
  for (int m=0; m<nargs; m++) {
    if (inds[m]>=0) maxbytes += args[m].dat->size;
  }

  // set blocksize and number of blocks;
  // adaptive size based on 48kB of shared memory

  int bsize   = part_size;   // blocksize
  if (bsize==0) bsize = (48*1024/(64*maxbytes))*64;
  int nblocks = (set->size-1)/bsize + 1;

  // enlarge OP_plans array if needed

  if (ip==OP_plan_max) {
    // printf("allocating more memory for OP_plans \n");
    OP_plan_max += 10;
    OP_plans = (op_plan *) realloc(OP_plans,OP_plan_max*sizeof(op_plan));
    if (OP_plans==NULL) {
      printf(" op_plan error -- error reallocating memory for OP_plans\n");
      exit(-1);  
    }
  }

  // allocate memory for new execution plan and store input arguments

  OP_plans[ip].dats      = (op_dat *)malloc(nargs*sizeof(op_dat));
  OP_plans[ip].idxs      = (int *)malloc(nargs*sizeof(int));
  OP_plans[ip].maps      = (op_map *)malloc(nargs*sizeof(op_map));
  OP_plans[ip].accs      = (op_access *)malloc(nargs*sizeof(op_access));

  OP_plans[ip].nthrcol   = (int *)malloc(nblocks*sizeof(int));
  OP_plans[ip].thrcol    = (int *)malloc(set->size*sizeof(int));
  OP_plans[ip].offset    = (int *)malloc(nblocks*sizeof(int));
  OP_plans[ip].ind_maps  = (int **)malloc(ninds*sizeof(int *));
  OP_plans[ip].ind_offs  = (int *)malloc(nblocks*ninds*sizeof(int));
  OP_plans[ip].ind_sizes = (int *)malloc(nblocks*ninds*sizeof(int));
  OP_plans[ip].nindirect = (int *)calloc(ninds,sizeof(int));
  OP_plans[ip].loc_maps  = (short **)malloc(nargs*sizeof(short *));
  OP_plans[ip].nelems    = (int *)malloc(nblocks*sizeof(int));
  OP_plans[ip].ncolblk   = (int *)calloc(set->size,sizeof(int)); // max possibly needed
  OP_plans[ip].blkmap    = (int *)calloc(nblocks,sizeof(int));

  for (int m=0; m<ninds; m++) {
    int count = 0;
    for (int m2=0; m2<nargs; m2++)
      if (inds[m2]==m) count++;
    OP_plans[ip].ind_maps[m] = (int *)malloc(count*set->size*sizeof(int));
  }

  for (int m=0; m<nargs; m++) {
    if (inds[m]>=0)
      OP_plans[ip].loc_maps[m] = (short *)malloc(set->size*sizeof(short));
    else
      OP_plans[ip].loc_maps[m] = NULL;

    OP_plans[ip].dats[m] = args[m].dat;
    OP_plans[ip].idxs[m] = args[m].idx;
    OP_plans[ip].maps[m] = args[m].map;
    OP_plans[ip].accs[m] = args[m].acc;
  }

  OP_plans[ip].name      = name;
  OP_plans[ip].set       = set;
  OP_plans[ip].nargs     = nargs;
  OP_plans[ip].ninds     = ninds;
  OP_plans[ip].part_size = part_size;
  OP_plans[ip].nblocks   = nblocks;
  OP_plans[ip].count     = 1;
    
  OP_plan_index++;

  // define aliases

  op_dat    *dats = OP_plans[ip].dats;
  int       *idxs = OP_plans[ip].idxs;
  op_map    *maps = OP_plans[ip].maps;
  op_access *accs = OP_plans[ip].accs;

  int  *offset    = OP_plans[ip].offset;
  int  *nelems    = OP_plans[ip].nelems;
  int **ind_maps  = OP_plans[ip].ind_maps;
  int  *ind_offs  = OP_plans[ip].ind_offs;
  int  *ind_sizes = OP_plans[ip].ind_sizes;
  int  *nindirect = OP_plans[ip].nindirect;

  // allocate working arrays

  uint **work;
  work = (uint **)malloc(ninds*sizeof(uint *));

  for (int m=0; m<ninds; m++) {
    int m2 = 0;
    while(inds[m2]!=m) m2++;

    work[m] = (uint *)malloc((maps[m2]->to)->size*sizeof(uint));
  }

  int *work2;
  work2 = (int *)malloc(nargs*bsize*sizeof(int));  // max possibly needed

  // process set one block at a time

  double total_colors = 0;

  for (int b=0; b<nblocks; b++) {
    int bs = MIN(bsize, set->size - b*bsize);

    offset[b] = b*bsize;    // offset for block
    nelems[b] = bs;         // size of block

    // loop over indirection sets

    for (int m=0; m<ninds; m++) {

      // build the list of elements indirectly referenced in this block

      int ne = 0;  // number of elements
      for (int m2=0; m2<nargs; m2++) {
        if (inds[m2]==m) {
          for (int e=b*bsize; e<b*bsize+bs; e++)
            work2[ne++] = maps[m2]->map[idxs[m2]+e*maps[m2]->dim];
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

      /*
      if (OP_diags > 5) {
        printf(" indirection set %d: ",m);
        for (int e=0; e<ne; e++) printf(" %d",work2[e]);
        printf(" \n");
      }
      */


      // store mapping and renumbered mappings in execution plan

      for (int e=0; e<ne; e++) {
        ind_maps[m][nindirect[m]++] = work2[e];
        work[m][work2[e]] = e;   // inverse mapping
      }

      for (int m2=0; m2<nargs; m2++) {
        if (inds[m2]==m) {
          for (int e=b*bsize; e<b*bsize+bs; e++)
            OP_plans[ip].loc_maps[m2][e] =
                  work[m][maps[m2]->map[idxs[m2]+e*maps[m2]->dim]];
	}
      }

      if (b==0) {
        ind_offs[m+b*ninds]  = 0;
        ind_sizes[m+b*ninds] = nindirect[m];
      }
      else {
        ind_offs[m+b*ninds]  = ind_offs[m+(b-1)*ninds]
                             + ind_sizes[m+(b-1)*ninds];
        ind_sizes[m+b*ninds] = nindirect[m]
                             - ind_offs[m+b*ninds];
      }
    }

    // now colour main set elements

    for (int e=b*bsize; e<b*bsize+bs; e++) OP_plans[ip].thrcol[e]=-1;

    int repeat  = 1;
    int ncolor  = 0;
    int ncolors = 0;

    while (repeat) {
      repeat = 0;

      for (int m=0; m<nargs; m++) {
        if (inds[m]>=0)
          for (int e=b*bsize; e<b*bsize+bs; e++)
            work[inds[m]][maps[m]->map[idxs[m]+e*maps[m]->dim]] = 0;  // zero out color array
      }

      for (int e=b*bsize; e<b*bsize+bs; e++) {
        if (OP_plans[ip].thrcol[e] == -1) {
          int mask = 0;
          for (int m=0; m<nargs; m++)
            if (inds[m]>=0 && accs[m]==OP_INC)
              mask |= work[inds[m]][maps[m]->map[idxs[m]+e*maps[m]->dim]]; // set bits of mask

          int color = ffs(~mask) - 1;   // find first bit not set
          if (color==-1) {              // run out of colors on this pass
            repeat = 1;
          }
          else {
            OP_plans[ip].thrcol[e] = ncolor+color;
            mask    = 1 << color;
            ncolors = MAX(ncolors, ncolor+color+1);

            for (int m=0; m<nargs; m++)
              if (inds[m]>=0 && accs[m]==OP_INC)
                work[inds[m]][maps[m]->map[idxs[m]+e*maps[m]->dim]] |= mask; // set color bit
          }
        }
      }

      ncolor += 32;   // increment base level
    }

    OP_plans[ip].nthrcol[b] = ncolors;  // number of thread colors in this block
    total_colors += ncolors;

    // if(ncolors>1) printf(" number of colors in this block = %d \n",ncolors);

    // reorder elements by color?

  }


  // colour the blocks, after initialising colors to 0

  int *blk_col;
  blk_col = (int *) malloc(nblocks*sizeof(int));
  for (int b=0; b<nblocks; b++) blk_col[b] = -1;

  int repeat  = 1;
  int ncolor  = 0;
  int ncolors = 0;

  while (repeat) {
    repeat = 0;

    for (int m=0; m<nargs; m++) {
      if (inds[m]>=0) 
        for (int e=0; e<(maps[m]->to)->size; e++)
          work[inds[m]][e] = 0;               // zero out color arrays
    }

    for (int b=0; b<nblocks; b++) {
      if (blk_col[b] == -1) {          // color not yet assigned to block
        int  bs   = MIN(bsize, set->size - b*bsize);
        uint mask = 0;

        for (int m=0; m<nargs; m++) {
          if (inds[m]>=0 && accs[m]==OP_INC) 
            for (int e=b*bsize; e<b*bsize+bs; e++)
              mask |= work[inds[m]][maps[m]->map[idxs[m]+e*maps[m]->dim]]; // set bits of mask
        }

        int color = ffs(~mask) - 1;   // find first bit not set
        if (color==-1) {              // run out of colors on this pass
          repeat = 1;
        }
        else {
          blk_col[b] = ncolor + color;
          mask    = 1 << color;
          ncolors = MAX(ncolors, ncolor+color+1);

          for (int m=0; m<nargs; m++) {
            if (inds[m]>=0 && accs[m]==OP_INC)
              for (int e=b*bsize; e<b*bsize+bs; e++)
                work[inds[m]][maps[m]->map[idxs[m]+e*maps[m]->dim]] |= mask;
          }
        }
      }
    }

    ncolor += 32;   // increment base level
  }

  // store block mapping and number of blocks per color


  OP_plans[ip].ncolors = ncolors;

  for (int b=0; b<nblocks; b++)
    OP_plans[ip].ncolblk[blk_col[b]]++;  // number of blocks of each color

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

  // reorder blocks by color?


  // work out shared memory requirements

  OP_plans[ip].nshared = 0;
  double total_shared = 0;

  for (int b=0; b<nblocks; b++) {
    int nbytes = 0;
    for (int m=0; m<ninds; m++) {
      int m2 = 0;
      while(inds[m2]!=m) m2++;

      nbytes += ROUND_UP(ind_sizes[m+b*ninds]*dats[m2]->size);
    }
    OP_plans[ip].nshared = MAX(OP_plans[ip].nshared,nbytes);
    total_shared += nbytes;
  }

  // work out total bandwidth requirements

  OP_plans[ip].transfer  = 0;
  OP_plans[ip].transfer2 = 0;
  double transfer3 = 0;

  for (int b=0; b<nblocks; b++) {
    for (int m=0; m<nargs; m++) {
      if (inds[m]<0) {
        double fac = 2.0f;
        if (accs[m]==OP_READ) fac = 1.0f;
        if (dats[m]!=NULL) {
        OP_plans[ip].transfer  += fac*nelems[b]*dats[m]->size;
        OP_plans[ip].transfer2 += fac*nelems[b]*dats[m]->size;
        transfer3              += fac*nelems[b]*dats[m]->size; 
	}
      }
      else {
        OP_plans[ip].transfer  += nelems[b]*sizeof(short);
        OP_plans[ip].transfer2 += nelems[b]*sizeof(short);
        transfer3              += nelems[b]*sizeof(short);
      }
    }
    for (int m=0; m<ninds; m++) {
      int m2 = 0;
      while(inds[m2]!=m) m2++;
      double fac = 2.0f;
      if (accs[m2]==OP_READ) fac = 1.0f;
      OP_plans[ip].transfer +=
                      fac*ind_sizes[m+b*ninds]*dats[m2]->size;

      // work out how many cache lines are used by indirect addressing

      int i_map, l_new, l_old;
      int e0 = ind_offs[m+b*ninds];
      int e1 = e0 + ind_sizes[m+b*ninds];

      l_old = -1;

      for (int e=e0; e<e1; e++) {
        i_map = ind_maps[m][e];
        l_new = (i_map*dats[m2]->size)/OP_cache_line_size;
        if (l_new>l_old) OP_plans[ip].transfer2 += fac*OP_cache_line_size;
        l_old = l_new;
        l_new = ((i_map+1)*dats[m2]->size-1)/OP_cache_line_size;
        OP_plans[ip].transfer2 += fac*(l_new-l_old)*OP_cache_line_size;
        l_old = l_new;
      }

      l_old = -1;

      for (int e=e0; e<e1; e++) {
        i_map = ind_maps[m][e];
        l_new = (i_map*dats[m2]->size)/(dats[m2]->dim*OP_cache_line_size);
        if (l_new>l_old) transfer3 += fac*dats[m2]->dim*OP_cache_line_size;
        l_old = l_new;
        l_new = ((i_map+1)*dats[m2]->size-1)/(dats[m2]->dim*OP_cache_line_size);
        transfer3 += fac*(l_new-l_old)*dats[m2]->dim*OP_cache_line_size;
        l_old = l_new;
      }

      // also include mappings to load/store data

      fac = 1.0f;
      if (accs[m2]==OP_RW) fac = 2.0f;
      OP_plans[ip].transfer  += fac*ind_sizes[m+b*ninds]*sizeof(int);
      OP_plans[ip].transfer2 += fac*ind_sizes[m+b*ninds]*sizeof(int);
      transfer3              += fac*ind_sizes[m+b*ninds]*sizeof(int);
    }
  }

  // print out useful information

  if (OP_diags>1) {
    printf(" number of blocks       = %d \n",nblocks);
    printf(" number of block colors = %d \n",OP_plans[ip].ncolors);
    printf(" maximum block size     = %d \n",bsize);
    printf(" average thread colors  = %.2f \n",total_colors/nblocks);
    printf(" shared memory required = %.2f KB \n",OP_plans[ip].nshared/1024.0f);
    printf(" average data reuse     = %.2f \n",maxbytes*(set->size/total_shared));
    printf(" data transfer (used)   = %.2f MB \n",
                                       OP_plans[ip].transfer/(1024.0f*1024.0f));
    printf(" data transfer (total)  = %.2f MB \n",
                                       OP_plans[ip].transfer2/(1024.0f*1024.0f));
    printf(" SoA/AoS transfer ratio = %.2f \n\n",
                                       transfer3 / OP_plans[ip].transfer2);
  }

  // validate plan info

  op_plan_check(OP_plans[ip],ninds,inds);

  // free work arrays

  for (int m=0; m<ninds; m++) free(work[m]);
  free(work);
  free(work2);
  free(blk_col);

  // return pointer to plan

  return &(OP_plans[ip]);
}


//
//include mpi halo functions
//
#include "op_mpi_core.c"


//
//include mpi partitioning functions
//
#include "op_mpi_part_core.c"
