
//
// headers
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "op_datatypes.h"


//
// op routine declarations
//

void op_init(int, char **);

void op_decl_set(int, int, float *, op_set &, char const *, op_ptr &, char const *);

void op_decl_ptr(op_set, op_set, int, int *, op_ptr &, char const *);

void op_decl_ddat(op_set, int, op_datatype, double *, op_dat &, char const *);

void op_decl_fdat(op_set, int, op_datatype, float *, op_dat &, char const *);

void op_decl_idat(op_set, int, op_datatype, int *, op_dat &, char const *);

void op_decl_dat(op_set set, int dim, op_datatype type, double *dat, op_dat &data, char const *name){
    op_decl_ddat(       set,     dim,             type,         dat,         data,             name);
}

void op_decl_dat(op_set set, int dim, op_datatype type, float *dat, op_dat &data, char const *name){
    op_decl_fdat(       set,     dim,             type,        dat,         data,             name);
}

void op_decl_dat(op_set set, int dim, op_datatype type, int *dat, op_dat &data, char const *name){
    op_decl_idat(       set,     dim,             type,      dat,         data,             name);
}

void op_fetch_data(op_dat);

void op_diagnostic_output();


//
// op_par_loop routine for 3 arguments -- need to reproduce for other numbers of arguments
//

template <class T0, class T1, class T2>
void op_par_loop_3(void (*kernel)(T0 *, T1 *, T2 *), char const * name, op_set set,
    op_dat arg0, int idx0, op_ptr ptr0, int dim0, op_datatype typ0, op_access acc0,
    op_dat arg1, int idx1, op_ptr ptr1, int dim1, op_datatype typ1, op_access acc1,
    op_dat arg2, int idx2, op_ptr ptr2, int dim2, op_datatype typ2, op_access acc2){

  int      nargs = 3;

  int         idxs[3] = {idx0, idx1, idx2};
  int         dims[3] = {dim0, dim1, dim2};
  op_datatype typs[3] = {typ0, typ1, typ2};
  op_access   accs[3] = {acc0, acc1, acc2};

  int    *ptr_ptrs[3] = {ptr0.ptr, ptr1.ptr, ptr2.ptr};
  int     ptr_dims[3] = {ptr0.dim, ptr1.dim, ptr2.dim};

  float    *f_dats[3] = {arg0.fdat, arg1.fdat, arg2.fdat};
  double   *d_dats[3] = {arg0.ddat, arg1.ddat, arg2.ddat};
  int      *i_dats[3] = {arg0.idat, arg1.idat, arg2.idat};

  float    *f_args[3];
  double   *d_args[3];
  int      *i_args[3];
  void       *args[3];

  int n2, ninds;

  // consistency checks

  if (OP_DIAGS>0) {
    int          ptr_from[3] = {ptr0.from.index, ptr1.from.index, ptr2.from.index};
    int            ptr_to[3] = {ptr0.to.index,   ptr1.to.index,   ptr2.to.index  };
    int          arg_sets[3] = {arg0.set.index,  arg1.set.index,  arg2.set.index };
    int          arg_dims[3] = {arg0.dim,        arg1.dim,        arg2.dim};
    op_datatype arg_types[3] = {arg0.type,       arg1.type,       arg2.type};

    for (int m=0; m<nargs; m++) {
      if (idxs[m]>=0) ninds++;

      if (set.index != ptr_from[m] || arg_sets[m] != ptr_to[m]) {
        printf("error: wrong pointer for arg %d in kernel \"%s\"\n",m,name); exit(1);
      }
      if (arg_types[m] != typs[m]) {
        printf("error: wrong datatype for arg %d in kernel \"%s\"\n",m,name); exit(1);
      }
      if (arg_dims[m] != dims[m]) {
        printf("error: wrong dimension for arg %d in kernel \"%s\"\n",m,name); exit(1);
      }
      if (ptr_dims[m] <= idxs[m]) {
        printf(" %d %d",ptr_dims[m],idxs[m]);
        printf("error: invalid pointer index for arg %d in kernel \"%s\"\n",m,name); exit(1);
      }
    }
  }

  if (OP_DIAGS>1) {
    if (ninds==0)
      printf(" kernel routine w/o indirection:  %s \n",name);
    else
      printf(" kernel routine with indirection: %s \n",name);
  }

  // allocate memory for local arrays, and
  // initialise to zero (important for increments)

  for (int m=0; m<nargs; m++) {
    if (typs[m] == OP_FLOAT) {
      args[m]   = calloc(dims[m],sizeof(float));
      f_args[m] = (float *) args[m];
    }
    else if (typs[m] == OP_DOUBLE) {
      args[m]   = calloc(dims[m],sizeof(double));
      d_args[m] = (double *) args[m];
    }
    else if (typs[m] == OP_INT) {
      args[m]   = calloc(dims[m],sizeof(int));
      i_args[m] = (int *) args[m];
    }
  }


  // loop over set elements

  for (int n=0; n<set.size; n++) {
    for (int m=0; m<nargs; m++) {
      if (accs[m]==OP_READ || accs[m]==OP_RW) {

        // if(OP_DIAGS>1 && n==0) printf(" m=%d, READ/RW read \n",m);

        if (ptr_dims[m]==0)                 // identity mapping
          n2 = n;
        else                                // standard pointers
          n2 = ptr_ptrs[m][idxs[m]+n*ptr_dims[m]];

        if (typs[m] == OP_FLOAT)
          for(int p=0; p<dims[m]; p++) f_args[m][p] = f_dats[m][p+n2*dims[m]];
        else if (typs[m] == OP_DOUBLE)
          for(int p=0; p<dims[m]; p++) d_args[m][p] = d_dats[m][p+n2*dims[m]];
        else if (typs[m] == OP_INT)
          for(int p=0; p<dims[m]; p++) i_args[m][p] = i_dats[m][p+n2*dims[m]];
      }
    }

    // call kernel function, passing in pointers to data
    
    kernel((T0 *)args[0],(T1 *)args[1],(T2 *)args[2]);

    for (int m=0; m<nargs; m++) {
      if (accs[m]==OP_WRITE || accs[m]==OP_RW) {

        // if(OP_DIAGS>1 && n==0) printf(" m=%d, WRITE/RW write \n",m);

        if (ptr_dims[m]==0)                 // identity mapping
          n2 = n;
        else                                // standard pointers
          n2 = ptr_ptrs[m][idxs[m]+n*ptr_dims[m]];

        if (typs[m] == OP_FLOAT)
          for(int p=0; p<dims[m]; p++) f_dats[m][p+n2*dims[m]] = f_args[m][p];
	else if (typs[m] == OP_DOUBLE)
          for(int p=0; p<dims[m]; p++) d_dats[m][p+n2*dims[m]] = d_args[m][p];
	else if (typs[m] == OP_INT)
          for(int p=0; p<dims[m]; p++) i_dats[m][p+n2*dims[m]] = i_args[m][p];
      }

      if (accs[m]==OP_INC) {

        // if(OP_DIAGS>1 && n==0) printf(" m=%d, INC increment \n",m);

        if (ptr_dims[m]==0)                 // identity mapping
          n2 = n;
        else                                // standard pointers
          n2 = ptr_ptrs[m][idxs[m]+n*ptr_dims[m]];

        if (typs[m] == OP_FLOAT)
          for(int p=0; p<dims[m]; p++) {
            f_dats[m][p+n2*dims[m]] += f_args[m][p];
            f_args[m][p] = 0.0f;
	  }
	else if (typs[m] == OP_DOUBLE)
          for(int p=0; p<dims[m]; p++) {
            d_dats[m][p+n2*dims[m]] += d_args[m][p];
            d_args[m][p] = 0.0;
	  }
        else if (typs[m] == OP_INT)
          for(int p=0; p<dims[m]; p++) {
            i_dats[m][p+n2*dims[m]] += i_args[m][p];
            i_args[m][p] = 0;
	  }
      }
    }
  }
}

