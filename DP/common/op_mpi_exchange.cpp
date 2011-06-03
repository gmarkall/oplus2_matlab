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
 
 
/*---------------------------MPI Halo Exchange -------------------------------*/

void exchange_halo(op_set set, op_dat dat, op_access acc, int idx)
{
    if((idx != -1) && (acc == OP_READ || acc == OP_RW ) && 
    	(dirtybit[dat.index] == 1))
    {
    	//printf("Exchanging Halo of data array %10s\n",dat.name);
   	
    	set_import_list imp_exec_list = *OP_import_sets_list[dat.set.index];
    	nonexec_set_import_list imp_nonexec_list = *OP_import_nonexec_sets_list[dat.set.index];
    	
    	set_export_list exp_exec_list = *OP_export_sets_list[dat.set.index];
    	nonexec_set_export_list exp_nonexec_list = *OP_export_nonexec_sets_list[dat.set.index];
    	
    	//-------first exchange exec elements related to this data array--------
    
    	//sanity checks
    	if(compare_sets(imp_exec_list.set,dat.set)==0) 
    	    { printf("Error: Import list and set mismatch\n"); exit(2);}
    	if(compare_sets(exp_exec_list.set,dat.set)==0) 
    	    {printf("Error: Export list and set mismatch\n"); exit(2);}
        	
    	MPI_Request request_send[exp_exec_list.ranks_size];
    	char* sbuf[exp_exec_list.ranks_size];
    	//prepare execute set element data to be exported
    	for(int i=0; i<exp_exec_list.ranks_size; i++) {
    	    sbuf[i] = allocate(sbuf[i],exp_exec_list.sizes[i]*dat.size);
    	    for(int j = 0; j<exp_exec_list.sizes[i]; j++)
    	    {
    	    	int set_elem_index = exp_exec_list.exp_list[exp_exec_list.disps[i]+j];
    	    	memcpy(&sbuf[i][j*dat.size],(void *)&dat.dat[dat.size*(set_elem_index)],dat.size);
    	    }
    	    //printf("export from %d to %d data %10s, number of elements of size %d | sending:\n ",
      	    //  	      my_rank, exp_exec_list.ranks[i], dat.name,exp_exec_list.sizes[i]);
    	    MPI_Isend(sbuf[i],  dat.size*exp_exec_list.sizes[i],  
    	    	MPI_CHAR, exp_exec_list.ranks[i],
    	    	dat.index, OP_MPI_WORLD, &request_send[i]);
    	}
    	
    	int init = dat.set.size*dat.size;
    	
    	for(int i=0; i<imp_exec_list.ranks_size; i++) {
    	    //printf("import on to %d from %d data %10s, number of elements of size %d | recieving:\n ",
      	    //  	  my_rank, imp_exec_list.ranks[i], dat.name, imp_exec_list.sizes[i]);
    	    MPI_Recv(&(OP_dat_list[dat.index]->dat[init+imp_exec_list.disps[i]*dat.size]),
    	        dat.size*imp_exec_list.sizes[i], MPI_CHAR, imp_exec_list.ranks[i], dat.index,
    	        OP_MPI_WORLD, MPI_STATUSES_IGNORE);
    	}
    	MPI_Waitall(exp_exec_list.ranks_size,request_send, MPI_STATUSES_IGNORE );
        for(int i=0; i<exp_exec_list.ranks_size; i++) free(sbuf[i]);
            

    	//-----second exchange nonexec elements related to this data array------
    	//sanity checks
   	if(compare_sets(imp_nonexec_list.set,dat.set)==0) 
   	    { printf("Error: Non-Import list and set mismatch"); exit(2);}
    	if(compare_sets(exp_nonexec_list.set,dat.set)==0) 
    	    {printf("Error: Non-Export list and set mismatch"); exit(2);}

    	
    	MPI_Request request_nonexec_send[exp_nonexec_list.ranks_size];
    	char* sbuf_nonexec[exp_nonexec_list.ranks_size];
    	//prepare execute set element data to be exported
    	for(int i=0; i<exp_nonexec_list.ranks_size; i++) {
    	    sbuf_nonexec[i] = allocate(sbuf_nonexec[i],exp_nonexec_list.sizes[i]*dat.size);
    	    for(int j = 0; j<exp_nonexec_list.sizes[i]; j++)
    	    {
    	    	int set_elem_index = exp_nonexec_list.exp_list[exp_nonexec_list.disps[i]+j];
    	    	memcpy(&sbuf_nonexec[i][j*dat.size],(void *)&dat.dat[dat.size*(set_elem_index)],dat.size);
    	    }
    	    
    	    MPI_Isend(sbuf_nonexec[i],  dat.size*exp_nonexec_list.sizes[i],  
    	    	MPI_CHAR, exp_nonexec_list.ranks[i],
    	    	dat.index, OP_MPI_WORLD, &request_nonexec_send[i]);
    	}
    	
    	int nonexec_init = (dat.set.size+imp_exec_list.size)*dat.size;
    	
    	for(int i=0; i<imp_nonexec_list.ranks_size; i++) {
    	    MPI_Recv(&(OP_dat_list[dat.index]->dat[nonexec_init+imp_nonexec_list.disps[i]*dat.size]),
    	        dat.size*imp_nonexec_list.sizes[i], MPI_CHAR, imp_nonexec_list.ranks[i], dat.index,
    	        OP_MPI_WORLD, MPI_STATUSES_IGNORE);
    	}
    	MPI_Waitall(exp_nonexec_list.ranks_size,request_nonexec_send, MPI_STATUSES_IGNORE );
        for(int i=0; i<exp_nonexec_list.ranks_size; i++) free(sbuf_nonexec[i]);
    	    	
    	//clear dirty bit
    	dirtybit[dat.index] = 0;
    }
}


void exchange_halo(op_set set, void* value, op_access acc, int idx)
{
    //not a data array for halo exchange - do nothing
}

void set_dirtybit(op_dat dat, op_access acc)
{
    if(acc == OP_INC || acc == OP_WRITE || acc == OP_RW)
    	dirtybit[dat.index] = 1;
}

void set_dirtybit(void* value, op_access acc)
{
    //not a data array for halo exchange - do nothing
}





