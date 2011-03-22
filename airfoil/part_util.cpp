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


/* 
 * written by: Gihan R. Mudalige, 01-03-2011
 */
 
 

/*---------utility functions for trivial block partitioning of sets-----------*/

int compute_local_size (int global_size, int mpi_comm_size, int mpi_rank )
{
  	  int local_size = global_size/mpi_comm_size;
  	  int remainder = (int)fmod(global_size,mpi_comm_size);
  
  	  if (mpi_rank < remainder)
  	  {
  	  	  local_size = local_size + 1;
  	  
  	  }
  	  return local_size;
}

int get_partition(int value, int* part_range, int* local_index, int comm_size)
{
	for(int i = 0; i<comm_size; i++)
  	{
  		if (value >= part_range[2*i] & value <= part_range[2*i+1])
  		{
  		    	*local_index = value -  part_range[2*i];
  		    	return i;
  		}
	}
	return 0;
}

void scatter_double_array(double* g_array, double* l_array, int comm_size, int g_size, 
	int l_size, int elem_size)
{
  	  int* sendcnts = (int *) malloc(comm_size*sizeof(int));
  	  int* displs = (int *) malloc(comm_size*sizeof(int));
  	  int disp = 0;
  
  	  for(int i = 0; i<comm_size; i++)
  	  {
  	  	  sendcnts[i] =   elem_size*compute_local_size (g_size, comm_size, i);	
  	  }
  	  for(int i = 0; i<comm_size; i++)
  	  {	
  	  	  displs[i] =   disp;
  	  	  disp = disp + sendcnts[i];
  	  }
  
  	  MPI_Scatterv(g_array, sendcnts, displs, MPI_DOUBLE, l_array, 
  	  	  l_size*elem_size, MPI_DOUBLE, 0,  MPI_COMM_WORLD ); 
}

void scatter_int_array(int* g_array, int* l_array, int comm_size, int g_size, 
	int l_size, int elem_size)
{
  	  int* sendcnts = (int *) malloc(comm_size*sizeof(int));
  	  int* displs = (int *) malloc(comm_size*sizeof(int));
  	  int disp = 0;
  
  	  for(int i = 0; i<comm_size; i++)
  	  {
  	  	  sendcnts[i] =   elem_size*compute_local_size (g_size, comm_size, i);	
  	  }
  	  for(int i = 0; i<comm_size; i++)
  	  {	
  	  	  displs[i] =   disp;
  	  	  disp = disp + sendcnts[i];
  	  }
  
  	  MPI_Scatterv(g_array, sendcnts, displs, MPI_INT, l_array, 
  	  	  l_size*elem_size, MPI_INT, 0,  MPI_COMM_WORLD ); 
}

void get_part_range(int** part_range, int comm_size)
{
	for(int n=0; n<OP_set_index; n++) {
		
		op_set set=*OP_set_list[n];
		   
		part_range[n] = (int *) malloc(2*comm_size*sizeof(int)); 
		int disp = 0;
		
		//discover global size of set - inefficient to do this via an all_reduce
		int global_size;
		MPI_Allreduce(&set.size, &global_size,1, MPI_INT, MPI_SUM , MPI_COMM_WORLD );
		
		printf("%10s has local size: %10d | global size: %10d\n",set.name,set.size,global_size);
		for(int i = 0; i<comm_size; i++){
		      
		      part_range[n][2*i] = disp;
		      disp = disp + compute_local_size (global_size, comm_size, i) - 1;
		      part_range[n][2*i+1] = disp;
		      disp++;
		      printf("range of %10s in partition %d: %d-%d\n",set.name,i,
			      part_range[n][2*i], part_range[n][2*i+1]);
		} 
	
	}
		
}
