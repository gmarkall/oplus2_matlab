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
enum op_datatype { OP_FLOAT, OP_DOUBLE, OP_INT, OP_BOOL };

//
// run-time type-checking routines
//

inline int type_error(const double *,op_datatype type){return (type != OP_DOUBLE);}
inline int type_error(const float  *,op_datatype type){return (type != OP_FLOAT);}
inline int type_error(const int    *,op_datatype type){return (type != OP_INT);}
inline int type_error(const bool   *,op_datatype type){return (type != OP_BOOL);}

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
  op_datatype type;   // datatype
  char const *name;   // name of dataset
} op_dat;

// null set
#define OP_NULL_SET (op_set) {0,0,"null"}

// identity mapping
#define OP_ID (op_ptr) {OP_NULL_SET,OP_NULL_SET,0,-1,NULL,"id"}

// global identifier
#define OP_GBL (op_ptr) {OP_NULL_SET,OP_NULL_SET,0,-2,NULL,"gbl"}

typedef struct {
  // input arguments
  char const  *name;
  int          set_index, nargs;
  int         *arg_idxs, *idxs, *ptr_idxs, *dims;
  op_datatype *typs;
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

#endif

//
//  min / max definitions
//

#ifndef MIN
#define MIN(a,b) ((a<b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a>b) ? (a) : (b))
#endif


