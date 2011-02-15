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
//     Nonlinear airfoil lift calculation
//
//     Written by Mike Giles, 2010, based on FORTRAN code
//     by Devendra Ghate and Mike Giles, 2005
//

//
// standard headers
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// global constants

float gam, gm1, cfl, eps, mach, alpha;

//
// OP header file
//

#include "op_seq.h"

//
// kernel routines for parallel loops
//

#include "input.h"
#include "save_soln.h"
#include "adt_calc.h"
#include "res_calc.h"
#include "update.h"

// main program

int main(int argc, char **argv){

  int    *ecell, *boun, *edge, *cell;
  float  *x, *q, *qold, *adt, *res;

  int    nnode,ncell,nedge, niter;
  float  rms;

  op_set nodes, edges, cells;
  op_ptr pedge, pecell, pcell;
  op_dat p_x, p_q, p_qold, p_res, p_adt, p_boun;

  // read in grid and flow data

  
  input(nnode,ncell,nedge,x,q,cell,edge,ecell,boun);
  qold = (float *) malloc(4*ncell*sizeof(float));
  res  = (float *) malloc(4*ncell*sizeof(float));
  adt  = (float *) malloc(  ncell*sizeof(float));

  // initialise residual

  for (int n=0; n<4*ncell; n++) res[n]=0.0;

  // OP initialisation

  op_init(argc,argv,2);

  // declare sets, pointers, datasets and global constants

  op_decl_set(nnode, nodes,"nodes");
  op_decl_set(nedge, edges,"edges");
  op_decl_set(ncell, cells,"cells");

  op_decl_ptr(edges,nodes,2,edge, pedge, "pedge");
  op_decl_ptr(edges,cells,2,ecell,pecell,"pecell");
  op_decl_ptr(cells,nodes,4,cell, pcell, "pcell");

  op_decl_dat(edges,1,"int"   ,boun,p_boun,"p_boun");
  op_decl_dat(nodes,2,"float",x   ,p_x   ,"p_x");
  op_decl_dat(cells,4,"float",q   ,p_q   ,"p_q");
  op_decl_dat(cells,4,"float",qold,p_qold,"p_qold");
  op_decl_dat(cells,1,"float",adt ,p_adt ,"p_adt");
  op_decl_dat(cells,4,"float",res ,p_res ,"p_res");

  op_decl_const(1,"float",&gam  ,"gam");
  op_decl_const(1,"float",&gm1  ,"gm1");
  op_decl_const(1,"float",&cfl  ,"cfl");
  op_decl_const(1,"float",&eps  ,"eps");
  op_decl_const(1,"float",&mach ,"mach");
  op_decl_const(1,"float",&alpha,"alpha");

  op_diagnostic_output();

// main time-marching loop

  niter = 1000;

  for(int iter=1; iter<=niter; iter++) {

//  save old flow solution

    op_par_loop_2(save_soln,"save_soln", cells,
                  p_q,   -1,OP_ID, 4,"float",OP_READ,
                  p_qold,-1,OP_ID, 4,"float",OP_WRITE);

//  predictor/corrector update loop

    for(int k=0; k<2; k++) {

//    calculate area/timstep

      op_par_loop_6(adt_calc,"adt_calc",cells,
                    p_x,   0,pcell, 2,"float",OP_READ,
                    p_x,   1,pcell, 2,"float",OP_READ,
                    p_x,   2,pcell, 2,"float",OP_READ,
                    p_x,   3,pcell, 2,"float",OP_READ,
                    p_q,  -1,OP_ID, 4,"float",OP_READ,
                    p_adt,-1,OP_ID, 1,"float",OP_WRITE);

//    calculate flux residual

      op_par_loop_9(res_calc,"res_calc",edges,
                    p_x,    0,pedge, 2,"float",OP_READ,
                    p_x,    1,pedge, 2,"float",OP_READ,
                    p_q,    0,pecell,4,"float",OP_READ,
                    p_q,    1,pecell,4,"float",OP_READ,
                    p_adt,  0,pecell,1,"float",OP_READ,
                    p_adt,  1,pecell,1,"float",OP_READ,
                    p_res,  0,pecell,4,"float",OP_INC,
                    p_res,  1,pecell,4,"float",OP_INC,
                    p_boun,-1,OP_ID ,1,"int",   OP_READ);

      //      printf("res[0-3] = %14.10e %14.10e %14.10e %14.10e \n",res[0],res[1],res[2],res[3]);

//    update flow field

      rms = 0.0;

      op_par_loop_5(update,"update",cells,
                    p_qold,-1,OP_ID, 4,"float",OP_READ,
                    p_q,   -1,OP_ID, 4,"float",OP_WRITE,
                    p_res, -1,OP_ID, 4,"float",OP_RW,
                    p_adt, -1,OP_ID, 1,"float",OP_READ,
                    &rms,  -1,OP_GBL,1,"float",OP_INC);
    }

//  print iteration history

    rms = sqrt(rms/(float) ncell);

    if (iter%100 == 0)
      printf(" %d  %10.5e \n",iter,rms);
  }

  op_timing_output();
}

