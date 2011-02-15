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


//  #define op_par_loop_3(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)  op_par_loop_##a0(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)



#ifndef OP_DIAGS
#define OP_DIAGS 6
#endif


//
// OP datatypes
//

enum op_datatype { OP_FLOAT, OP_DOUBLE, OP_INT };
enum op_access   { OP_READ, OP_WRITE, OP_RW, OP_INC };

typedef struct {
  int         size, dim, index;
  float      *x;
  char const *name;
} op_set;


typedef struct {
  op_set      from, to;
  int         dim, index, *ptr; 
  char const *name;
} op_ptr;


typedef struct {
  op_set      set;
  int         dim, index;
  float      *fdat, *fdat_d;
  double     *ddat, *ddat_d;
  int        *idat, *idat_d;
  op_datatype type;
  char const *name;
} op_dat;


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


//
//  min / max definitions
//

#ifndef MIN
#define MIN(a,b) ((a<b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a>b) ? (a) : (b))
#endif


