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
 

 
/*-------------------------------Debug functions------------------------------*/

//print all lists for a given set 
void print_set_list(op_set set, int rank, int exec, int import)
{
    //exec import list
    if(exec==1 & import ==1)
    {
    	set_import_list l=*OP_import_sets_list[set.index];
    	printf("%d Exec import list for set %10s has %d neighbors\n",rank, set.name, l.ranks_size);
    	for(int j=0; j<l.ranks_size;j++)
    	{
    	    for(int k=0;k<l.sizes[j];k++)
    	    printf("On rank %d import from %d set %6s element %d\n",
    	    	rank,l.ranks[j],l.set.name,l.imp_list[l.disps[j]+k]);
    	}
    }
    //exec export list
    if(exec==1 & import ==0)
    {
    	set_export_list l=*OP_export_sets_list[set.index];
    	printf("On rank %d Exec export list for set %10s\n",rank, set.name);
    	for(int j=0; j<l.ranks_size;j++)
    	{
    	    for(int k=0;k<l.sizes[j];k++)
    	    printf("export to %d on rank %d set %6s contains %d\n",
    	    	l.ranks[j],rank,l.set.name,l.exp_list[l.disps[j]+k]);
    	}
    }
    
    //nonexec import list
    if(exec==0 & import ==1)
    {
    	nonexec_set_import_list l=*OP_import_nonexec_sets_list[set.index];
    	printf("On rank %d NonExec import list for set %10s\n",rank, set.name);
    	for(int j=0; j<l.ranks_size;j++)
    	{
    	    for(int k=0;k<l.sizes[j];k++)
    	    printf("On rank %d import from %d set %6s element %d\n",
    	    	rank,l.ranks[j],l.set.name,l.imp_list[l.disps[j]+k]);
    	}
    }
    
    //nonexec export list
    if(exec==0 & import ==0)
    {
    	nonexec_set_export_list l=*OP_export_nonexec_sets_list[set.index];
    	printf("On rank %d NonExec export list for set %10s\n",rank, set.name);
    	for(int j=0; j<l.ranks_size;j++)
    	{
    	    for(int k=0;k<l.sizes[j];k++)
    	    printf("export to %d on rank %d set %6s contains %d\n",
    	    	l.ranks[j],rank,l.set.name,l.exp_list[l.disps[j]+k]);
    	}
    }
}

//print a mapping table completely (either original or original+exec_list)
void print_map(op_map map, int rank, int length)
{
    printf("On rank %d Mapping table %10s\n",rank, map.name);
    
    for(int i = 0; i<length; i++) //for each element of this mapping table
    {
    	 printf("On rank %d map %10s contains|%d :",rank, map.name,i);
    	 for(int j=0; j<map.dim; j++) { //for each element pointed at by this entry
    	     
    	     printf("%d ", map.map[i*map.dim+j]);
    	 }
    	 printf("\n");    		
    }
}

//print a mapping table's import and/or export list
void print_map_list(op_map map, int rank, int import)
{
    //exec import list
    if(import ==1)
    {
    	map_import_list l=*OP_import_maps_list[map.index];
    	printf("On rank %d Exec import list for map %10s\n",rank, map.name);
    	for(int j=0; j<l.ranks_size;j++)
    	{
    	    for(int k=0;k<l.sizes[j];k++)
    	    printf("On rank %d import from %d map %6s entry %d\n",
    	    	rank, l.ranks[j], l.map.name, l.imp_list[l.disps[j]+k]);
    	}
    }
    
    //exec export list
    if(import ==0)
    {
    	map_export_list l=*OP_export_maps_list[map.index];
    	printf("On rank %d Exec export list for map %10s\n",rank, map.name);
    	for(int j=0; j<l.ranks_size;j++)
    	{
    	    for(int k=0;k<l.sizes[j];k++)
    	    printf("export to %d from rank %d  map %6s entry %d\n",
    	    	l.ranks[j], rank, l.map.name, l.exp_list[l.disps[j]+k]);
    	}
    }
}


//print a data array (i.e. proper element data values)
void print_dat_array(op_dat dat, int rank)
{
    op_dat data=*OP_dat_list[dat.index];
    printf("on rank %d: data array named %10s\n",
    	rank, dat.name);
    size_t mult = dat.size/dat.dim;
    
    for(int j = 0;j<dat.set.size;j++)
    {
    	double temp;
    	for(int e = 0; e<dat.dim;e++)
    	{
    	    memcpy (&temp, (void *)&(OP_dat_list[dat.index]->dat[(j*dat.dim+e)*mult]), mult);
    	    printf("%f ",temp);
    	}
    	printf("\n");
    }
}

//mpi_gather a data array (of type double) and print its values on proc 0
//(i.e. proper element data values)
void gatherprint_tofile(op_dat dat, int rank, int comm_size, int g_size, int l_size)
{
    op_dat data=*OP_dat_list[dat.index];
    
    if(l_size != data.set.size)
    {    printf("set size and data local size mismatch\n"); exit(1);}

    size_t mult = dat.size/dat.dim;
    
    double *l_array  = (double *) malloc(dat.dim*(dat.set.size)*sizeof(double));
    memcpy (l_array, (void *)&(OP_dat_list[dat.index]->dat[0]), dat.size*dat.set.size);
    
    
    int elem_size = dat.dim;
    int* recevcnts = (int *) malloc(comm_size*sizeof(int));
    int* displs = (int *) malloc(comm_size*sizeof(int));
    int disp = 0;
    double *g_array;
    
    if(rank==0) g_array  = (double *) malloc(elem_size*g_size*sizeof(double));
    
    for(int i = 0; i<comm_size; i++)
    {
    	  recevcnts[i] =   elem_size*compute_local_size (g_size, comm_size, i);	
    }
    for(int i = 0; i<comm_size; i++)
    {
    	displs[i] =   disp;
    	disp = disp + recevcnts[i];
    }
    
    MPI_Gatherv(l_array, l_size*elem_size, MPI_DOUBLE, g_array, recevcnts, 
  	      displs, MPI_DOUBLE, 0, OP_MPI_WORLD);
    
    
    if(rank==0)
    {
    	FILE *fp;
    	if ( (fp = fopen("out_grid.dat","w")) == NULL) {
    	    printf("can't open file out_grid.dat\n"); exit(-1);
    	}
    	
    	if (fprintf(fp,"%d %d\n",g_size, elem_size)<0)
    	{
    	    printf("error writing to out_grid.dat\n"); exit(-1);
    	}
    	
    	for(int i = 0; i< g_size; i++)
    	{
    	    for(int j = 0; j < elem_size; j++ )
    	    {
    	    	if (fprintf(fp,"%lf ",g_array[i*elem_size+j])<0)
    	    	{
    	    	    printf("error writing to out_grid.dat\n"); exit(-1);
    	    	}
    	    }
    	    fprintf(fp,"\n");
    	}
    	fclose(fp);
    }
}


//mpi_gather a data array (of type double) and print its values on proc 0 *TO BINARY FILE*
//(i.e. proper element data values)
void gatherprint_bin_tofile(op_dat dat, int rank, int comm_size, int g_size, int l_size)
{
    op_dat data=*OP_dat_list[dat.index];
    
    if(l_size != data.set.size)
    {    printf("set size and data local size mismatch\n"); exit(1);}

    size_t mult = dat.size/dat.dim;
    
    double *l_array  = (double *) malloc(dat.dim*(dat.set.size)*sizeof(double));
    memcpy (l_array, (void *)&(OP_dat_list[dat.index]->dat[0]), dat.size*dat.set.size);
    
    
    int elem_size = dat.dim;
    int* recevcnts = (int *) malloc(comm_size*sizeof(int));
    int* displs = (int *) malloc(comm_size*sizeof(int));
    int disp = 0;
    double *g_array;
    
    if(rank==0) g_array  = (double *) malloc(elem_size*g_size*sizeof(double));
    
    for(int i = 0; i<comm_size; i++)
    {
    	  recevcnts[i] =   elem_size*compute_local_size (g_size, comm_size, i);	
    }
    for(int i = 0; i<comm_size; i++)
    {
    	displs[i] =   disp;
    	disp = disp + recevcnts[i];
    }
    
    MPI_Gatherv(l_array, l_size*elem_size, MPI_DOUBLE, g_array, recevcnts, 
  	      displs, MPI_DOUBLE, 0, OP_MPI_WORLD);
    
    
    if(rank==0)
    {
    	FILE *fp;
    	if ( (fp = fopen("out_grid.bin","wb")) == NULL) {
    	    printf("can't open file out_grid.bin\n"); exit(-1);
    	}
    	
    	if (fwrite(&g_size, sizeof(int),1, fp)<0)
    	{
    	    printf("error writing to out_grid.bin"); exit(-1);
    	}
    	if (fwrite(&elem_size, sizeof(int),1, fp)<0)
    	{
    	    printf("error writing to out_grid.bin\n"); exit(-1);
    	}
    	
    	for(int i = 0; i< g_size; i++)
    	{
    	    if (fwrite( &g_array[i*elem_size], sizeof(double), elem_size, fp ) != 4)
    	    {
    	    	printf("error writing to out_grid.bin\n"); exit(-1);
    	    }
    	}
    	fclose(fp);
    }
}


void check_renumbering(int my_rank)
{
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=*OP_set_list[s];
    	
    	for(int m=0; m<OP_map_index; m++) { //for each maping table
    	    op_map map=*OP_map_list[m];
    	    
    	    if(compare_sets(map.to,set)==1) //need to select mappings TO this set
    	    {
    	    	map_import_list exec_map_list=*OP_import_maps_list[map.index];
    	    	set_import_list exec_set_list=*OP_import_sets_list[set.index];
    	    	nonexec_set_import_list nonexec_set_list=*OP_import_nonexec_sets_list[set.index];
    	    	
    	    	//for each entry in this mapping table: original+execlist
    	    	int len = map.from.size+exec_map_list.size;
    	    	for(int e = 0; e<len; e++)
    	    	{
    	    	    int part;
    	    	    int local_index;
    	    	    for(int j=0; j<map.dim; j++) { //for each element pointed at by this entry
    	    	    	if(map.map[e*map.dim+j]>= map.to.size +
    	    	    	    exec_set_list.size + nonexec_set_list.size)
    	    	    	printf("ERROR: Set %10s eleemnt %d out of bounds on %d \n",
    	    	    	    set.name, map.map[e*map.dim+j], my_rank);
    	    	    }
    	    	}
    	    }
    	}
    }
}


