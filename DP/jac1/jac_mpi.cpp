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
// test program for new OPlus2 development
//

//
// standard headers
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//mpi header
#include <mpi.h>

// global constants

double alpha;

//
// OP header file
//
#include "op_datatypes.h"
#include "op_mpi_seq.cpp" 

//
//trivial partitioning function
//
#include "part_util.cpp" 



//
//OP mpi halo creation and halo exchange functions
//
#include "op_mpi_list.cpp"
#include "op_mpi_exchange.cpp"
#include "op_mpi_seq.h" 



//
// kernel routines for parallel loops
//

#include "res.h"
#include "update.h"


// define problem size

#define NN       6
#define NITER    2


// main program

int main(int argc, char **argv){

  int my_rank;
  int comm_size;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  
  //timer
  double cpu_t1, cpu_t2, wall_t1, wall_t2;                                        
  double time;
  double max_time;
  
  int *pp;
  double *A, *r, *u, *du;
  
  int   nnode, nedge, n, e;
  double dx;

  op_set nodes, edges;
  op_map ppedge;
  op_dat p_A, p_r, p_u, p_du;

  /**------------------------BEGIN I/O and PARTITIONING ---------------------**/
  
  int g_nnode, g_nedge, g_dx, g_n, g_e;
  
  g_nnode = (NN-1)*(NN-1);
  g_nedge = (NN-1)*(NN-1) + 4*(NN-1)*(NN-2);
  g_dx    = 1.0f / ((double) NN);

  int *g_pp;
  double *g_A, *g_r, *g_u, *g_du;
  
  if(my_rank == 0)
  {
      printf("Global number of nodes, edges = %d, %d\n",g_nnode,g_nedge);
      
      g_pp = (int *)malloc(sizeof(int)*2*g_nedge);
    
      g_A  = (double *)malloc(sizeof(double)*g_nedge);
      g_r  = (double *)malloc(sizeof(double)*g_nnode);
      g_u  = (double *)malloc(sizeof(double)*g_nnode);
      g_du = (double *)malloc(sizeof(double)*g_nnode);
      
      // create matrix and r.h.s., and set coordinates needed for renumbering / partitioning
    
      g_e = 0;
    
      for (int i=1; i<NN; i++) {
	for (int j=1; j<NN; j++) {
	  g_n         = i-1 + (j-1)*(NN-1);
	  g_r[g_n]      = 0.0f;
	  g_u[g_n]      = 0.0f;
	  g_du[g_n]     = 0.0f;
    
	  g_pp[2*g_e]   = g_n;
	  g_pp[2*g_e+1] = g_n;
	  g_A[g_e]      = -1.0f;
	  g_e++;
    
	  for (int pass=0; pass<4; pass++) {
	    int i2 = i;
	    int j2 = j;
	    if (pass==0) i2 += -1;
	    if (pass==1) i2 +=  1;
	    if (pass==2) j2 += -1;
	    if (pass==3) j2 +=  1;
    
	    if ( (i2==0) || (i2==NN) || (j2==0) || (j2==NN) ) {
	      g_r[g_n] += 0.25f;
	    }
	    else {
	      g_pp[2*g_e]   = g_n;
	      g_pp[2*g_e+1] = i2-1 + (j2-1)*(NN-1);
	      g_A[g_e]      = 0.25f;
	      g_e++;
	    }
	  }
	}
      }
  }
  
  /* Compute local sizes */ 
  nnode = compute_local_size (g_nnode, comm_size, my_rank);
  nedge = compute_local_size (g_nedge, comm_size, my_rank);
  printf("Number of nodes, edges on process %d = %d, %d\n"
  	  ,my_rank,nnode,nedge);
  
  /*Allocate memory to hold local sets, mapping tables and data*/
  pp = (int *)malloc(2*sizeof(int)*nedge);
  
  A      = (double *) malloc(nedge*sizeof(double));
  r      = (double *) malloc(nnode*sizeof(double));
  u      = (double *) malloc(nnode*sizeof(double));
  du      = (double *) malloc(nnode*sizeof(double));
  
  /* scatter sets, mappings and data on sets*/
  scatter_int_array(g_pp, pp, comm_size, g_nedge,nedge, 2);
  scatter_double_array(g_A, A, comm_size, g_nedge,nedge, 1);
  scatter_double_array(g_r, r, comm_size, g_nnode,nnode, 1);
  scatter_double_array(g_u, u, comm_size, g_nnode,nnode, 1);
  scatter_double_array(g_du, du, comm_size, g_nnode,nnode, 1);
      
  if(my_rank == 0)
  { 	  /*Freeing memory allocated to gloabal arrays on rank 0 
  	  after scattering to all processes*/
	  free(g_pp);
	  free(g_A);
	  free(g_r);
	  free(g_u);
	  free(g_du);
  }
  
  /**------------------------END I/O and PARTITIONING ---------------------**/

  // OP initialisation

  op_init(argc,argv,2);

  // declare sets, pointers, and datasets

  op_decl_set(nnode, nodes,"nodes");
  op_decl_set(nedge, edges,"edges");

  op_decl_map(edges,nodes,2,pp, ppedge,"ppedge");

  op_decl_dat(edges,1,"double", A,  p_A, "p_A" );
  op_decl_dat(nodes,1,"double", r,  p_r, "p_r" );
  op_decl_dat(nodes,1,"double", u,  p_u, "p_u" );
  op_decl_dat(nodes,1,"double", du, p_du,"p_du");

  alpha = 1.0f;
  op_decl_const(1,"double",&alpha,"alpha");

  op_diagnostic_output();
  
  //mpi import/export list creation: should be included by the auto-generator
  op_list_create();
  
  //initialise timers for total execution wall time                                                         
  timers(&cpu_t1, &wall_t1); 
  

  // main iteration loop

  double u_sum, u_max, beta = 1.0f;

  for (int iter=0; iter<NITER; iter++) {
    op_par_loop_4(res,"res", edges,
                  p_A,  -1,OP_ID,  1,"double", OP_READ,
                  p_u,   1,ppedge, 1,"double", OP_READ,
                  p_du,  0,ppedge, 1,"double", OP_INC,
                  &beta,-1,OP_GBL, 1,"double", OP_READ);
    
    u_sum = 0.0f;
    u_max = 0.0f;
    op_par_loop_5(update,"update", nodes,
                  p_r,   -1,OP_ID, 1,"double",OP_READ,
                  p_du,  -1,OP_ID, 1,"double",OP_RW,
                  p_u,   -1,OP_ID, 1,"double",OP_INC,
                  &u_sum,-1,OP_GBL,1,"double",OP_INC,
                  &u_max,-1,OP_GBL,1,"double",OP_MAX);
    
    if(my_rank == 0)
    	printf("\n u max/rms = %f %f \n\n",u_max, sqrt(u_sum/g_nnode));
  }
  
  timers(&cpu_t2, &wall_t2); 
  
  //gatherprint_tofile(p_u, my_rank, comm_size, g_nnode, nnode);
  
  
  // print out results

  //printf("\n  Results after %d iterations:\n\n",NITER);

  //op_timing_output(my_rank);

  //print total time for niter interations 
  time = wall_t2-wall_t1;
  MPI_Reduce(&time,&max_time,1,MPI_DOUBLE, MPI_MAX,0, MPI_COMM_WORLD);
  if(my_rank==0)printf("Max total runtime = %f\n",max_time);
  
  MPI_Finalize();
  
  op_exit();
}
