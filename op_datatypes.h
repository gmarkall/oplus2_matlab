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


#ifndef OP_DIAGS
#define OP_DIAGS 6
#endif


//
// OP datatypes
//


#ifndef OP_DATATYPES
#define OP_DATATYPES

enum op_access   { OP_READ, OP_WRITE, OP_RW, OP_INC, OP_MIN, OP_MAX };

//
// run-time type-checking routines
//

typedef long long ll;
typedef unsigned long long ull;

inline int type_error(const double *,const char *type){return strcmp(type,"double");}
inline int type_error(const float  *,const char *type){return strcmp(type,"float" );}
inline int type_error(const int    *,const char *type){return strcmp(type,"int"   );}
inline int type_error(const uint   *,const char *type){return strcmp(type,"uint"  );}
inline int type_error(const ll     *,const char *type){return strcmp(type,"ll"    );}
inline int type_error(const ull    *,const char *type){return strcmp(type,"ull"   );}
inline int type_error(const bool   *,const char *type){return strcmp(type,"bool"  );}

//
// add in user's datatypes
//

#ifdef OP_USER_DATATYPES
#include <OP_USER_DATATYPES>
#endif

//
// zero constants
//

#define ZERO_double  0.0;
#define ZERO_float   0.0f;
#define ZERO_int     0;
#define ZERO_uint    0;
#define ZERO_ll      0;
#define ZERO_ull     0;
#define ZERO_bool    0;

//
// structures
//

typedef struct {
  int         size,   // number of elements in set
              index;  // index into list of sets 
  char const *name;   // name of set
} op_set;

typedef struct {
  op_set      from,   // set pointed from
              to;     // set pointed to
  int         dim,    // dimension of pointer
              index,  // index into list of pointers
             *ptr;    // array defining pointer
  char const *name;   // name of pointer
} op_ptr;

typedef struct {
  op_set      set;    // set on which data is defined
  int         dim,    // dimension of data
              index,  // index into list of datasets
              size;   // size of each element in dataset
  char       *dat,    // data on host
             *dat_d;  // data on device (GPU)
  char const *type,   // datatype
             *name;   // name of dataset
} op_dat;

// identity mapping
#define OP_ID  (op_ptr) {{0,0,"null"},{0,0,"null"},0,-1,NULL,"id"}

// global identifier
#define OP_GBL (op_ptr) {{0,0,"null"},{0,0,"null"},0,-2,NULL,"gbl"}

typedef struct {
  // input arguments
  char const  *name;
  int          set_index, nargs;
  int         *arg_idxs, *idxs, *ptr_idxs, *dims;
  char const **typs;
  op_access   *accs;

  // execution plan
  int        *nthrcol;  // number of thread colors for each block
  int        *thrcol;   // thread colors
  int        *offset;   // offset for primary set
  int       **ind_ptrs; // pointers for indirect datasets
  int       **ind_offs; // offsets for indirect datasets
  int       **ind_sizes;// sizes for indirect datasets
  int       **ptrs;     // regular pointers, renumbered as needed
  int        *nelems;   // number of elements in each block
  int         ncolors;  // number of block colors
  int        *ncolblk;  // number of blocks for each color
  int        *blkmap;   // block mapping
  int         nshared;  // bytes of shared memory required
} op_plan;


//
//  min / max definitions
//

#define MIN(a,b) ((a<b) ? (a) : (b))
#define MAX(a,b) ((a>b) ? (a) : (b))

//
// alignment macro based on example on page 50 of CUDA Programming Guide version 3.0
// rounds up to nearest multiple of 8 bytes
//

#define ROUND_UP(bytes) (((bytes) + 7) & ~7)

//
// OP function prototypes
//

void op_init(int, char **);

void op_decl_set(int, op_set &, char const *);

void op_decl_ptr(op_set, op_set, int, int *, op_ptr &, char const *);

void op_decl_dat_char(op_set, int, char const *, int, char *, op_dat &, char const *);

template < class T >
void op_decl_dat(op_set set, int dim, char const *type, T *dat, op_dat &data, char const *name){
  if (type_error(dat,type)) {
    printf("incorrect type specified for dataset \"%s\" \n",name);  exit(1);
  }
  op_decl_dat_char(set, dim, type, sizeof(T), (char *)dat, data, name);
}

void op_decl_const_char(int, char const *, int, char *, char const *);

template < class T >
void op_decl_const(int dim, char const *type, T *dat, char const *name){
  if (type_error(dat,type)) {
    printf("incorrect type specified for constant \"%s\" \n",name);  exit(1);
  }
  op_decl_const_char(dim, type, sizeof(T), (char *)dat, name);
}

void op_fetch_data(op_dat);

void op_diagnostic_output();

void op_exit();

#endif

