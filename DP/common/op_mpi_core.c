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

//mpi header
#include <mpi.h>

#include "op_mpi_core.h"

//utility functions
#include "mpi_util.c"


//
//MPI Communicator for halo creation and exchange
//
MPI_Comm OP_MPI_WORLD;


/**-----------------------MPI related global variables ----------------------**/
/*int OP_exp_maps_index=0;  //number of mapping table halo lists (export)
int OP_exp_exset_index=0; //number of halo lists with execute set elements (export) 

int OP_imp_maps_index=0; //number of mapping table halo lists (import)
int OP_imp_exset_index=0; //number halo lists with execute set elements (import) 

int OP_imp_nxset_index=0; //number halo lists with non-execute set elements (import) 
int OP_exp_nxset_index=0; //number halo lists with non-execute set elements (export) */



map_halo_list* OP_export_maps_list; 
set_halo_list* OP_export_sets_list; 

map_halo_list* OP_import_maps_list; 
set_halo_list* OP_import_sets_list;

set_halo_list* OP_import_nonexec_sets_list;
set_halo_list* OP_export_nonexec_sets_list;  

//global array to hold dirty_bits for op_dats
int* dirtybit;


/**-----------------------MPI related utility functions ---------------------**/

void get_part_range(int** part_range, int my_rank, int comm_size, MPI_Comm Comm)
{
    for(int s=0; s<OP_set_index; s++) {
    	op_set set=OP_set_list[s];
    	   	
    	int* sizes = (int *)malloc(sizeof(int)*comm_size);
    	if(part_range == NULL) {
    	    printf(" get_part_range -- error allocating memory: int* sizes\n");
    	    exit(-1);
    	}
    
    	MPI_Allgather(&set->size, 1, MPI_INT, sizes, 1, MPI_INT, Comm);
    	
    	part_range[set->index] = (int *) malloc(2*comm_size*sizeof(int));
    	if(part_range[set->index] == NULL) {
    	    printf(" get_part_range -- error allocating memory: part_range[s]\n");
    	    exit(-1);
    	}
    	
    	int disp = 0;
    	for(int i = 0; i<comm_size; i++){
    	    part_range[set->index][2*i] = disp;
    	    disp = disp + sizes[i] - 1;
    	    part_range[set->index][2*i+1] = disp;
    	    disp++;
    	    #if DEBUG
    	    if(my_rank == 0)
    	    printf("range of %10s in rank %d: %d-%d\n",set->name,i,
    	    	part_range[set->index][2*i], part_range[set->index][2*i+1]);
    	    #endif
    	}
    	free(sizes);
    }	
}

//get partition (i.e. mpi rank) where global_index is located and its local index
int get_partition(int global_index, int* part_range, int* local_index, int comm_size)
{
	for(int i = 0; i<comm_size; i++)
  	{
  		if (global_index >= part_range[2*i] & global_index <= part_range[2*i+1])
  		{
  		    	*local_index = global_index -  part_range[2*i];
  		    	return i;
  		}
	}
	return 0;
}

//convert a local index in to a global index
int get_global_index(int local_index, int partition, int* part_range, int comm_size)
{
    int g_index = part_range[2*partition]+local_index;
    if(g_index > part_range[2*(comm_size-1)+1]) printf("Global index larger than set size\n");
    return g_index;
}


void find_neighbors_map(map_halo_list List, int* neighbors, int* sizes, 
	int* ranks_size, int my_rank, int comm_size, MPI_Comm Comm)
{
    int *temp = malloc(comm_size*sizeof(int));
    int *r_temp = malloc(comm_size*comm_size*sizeof(int));
    
    for(int r = 0;r<comm_size*comm_size;r++)r_temp[r] = -99;
    for(int r = 0;r<comm_size;r++)temp[r] = -99;    
    
    int n = 0;
    
    for(int r =0; r<comm_size; r++)
    {
    	if(List->ranks[r]>=0) temp[List->ranks[r]] = List->sizes[r];
    }
    
    MPI_Allgather( temp, comm_size, MPI_INT, r_temp,
    	comm_size,MPI_INT,Comm);

    for(int i=0; i<comm_size; i++)
    {
    	if(i != my_rank)
    	{
    	    if( r_temp[i*comm_size+my_rank] > 0)
    	    {
    	    	neighbors[n] = i;
    	    	sizes[n] = r_temp[i*comm_size+my_rank];
    	    	n++;
    	    }    
    	}
    }
    *ranks_size = n;
    free(temp);free(r_temp);
}

void find_neighbors_set(set_halo_list List, int* neighbors, int* sizes, 
	int* ranks_size, int my_rank, int comm_size, MPI_Comm Comm)
{
    int *temp = malloc(comm_size*sizeof(int));
    int *r_temp = malloc(comm_size*comm_size*sizeof(int));
    
    for(int r = 0;r<comm_size*comm_size;r++)r_temp[r] = -99;
    for(int r = 0;r<comm_size;r++)temp[r] = -99;    
    
    int n = 0;
    
    for(int r =0; r<comm_size; r++)
    {
    	if(List->ranks[r]>=0) temp[List->ranks[r]] = List->sizes[r];
    }
    
    MPI_Allgather( temp, comm_size, MPI_INT, r_temp,
    	comm_size,MPI_INT,Comm);

    for(int i=0; i<comm_size; i++)
    {
    	if(i != my_rank)
    	{
    	    if( r_temp[i*comm_size+my_rank] > 0)
    	    {
    	    	neighbors[n] = i;
    	    	sizes[n] = r_temp[i*comm_size+my_rank];
    	    	n++;
    	    }    
    	}
    }
    *ranks_size = n;
    free(temp);free(r_temp);
}


void create_list(int* list, int* ranks, int* disps, int* sizes, int* ranks_size, int* total,
    int* temp_list, int size, int comm_size, int my_rank)
{
    int index = 0;
    int total_size = 0;
    
    //negative values set as an initialisation
    for(int r = 0;r<comm_size;r++)
    {
    	disps[r] = ranks[r] = -99;
    	sizes[r] = 0;
    }
    for(int r = 0;r<comm_size;r++)
    {
    	sizes[index] = disps[index] = 0;
    	
    	int* temp = (int *) malloc((size/2)*sizeof(int));
    	for(int i = 0;i<size;i=i+2)
    	{
    	    if(temp_list[i]==r)
    	    	temp[sizes[index]++] = temp_list[i+1];    	    	
    	}
    	if(sizes[index]>0)
    	{
    	    ranks[index] = r;
    	    //sort temp,
    	    quickSort(temp,0,sizes[index]-1);
    	    //eliminate duplicates in temp
    	    sizes[index] = removeDups(temp, sizes[index]);
    	    total_size = total_size + sizes[index];
    	    
    	    if(index > 0)
    	    	disps[index] = disps[index-1] +  sizes[index-1];
    	    //add to end of exp_list
    	    for(int e = 0;e<sizes[index];e++)
    	    	list[disps[index]+e] = temp[e];
    	    
    	    index++;
    	}
	free(temp);
    }
    
    *total = total_size;
    *ranks_size = index;
}

void create_set_export_list(op_set set, int* temp_list, int size, int comm_size, int my_rank)
{
    int* ranks = (int *) malloc(comm_size*sizeof(int));
    int* list = (int *) malloc((size/2)*sizeof(int));
    int* disps = (int *) malloc(comm_size*sizeof(int));
    int* sizes = (int *) malloc(comm_size*sizeof(int));
    
    int ranks_size = 0;
    int total_size = 0;  
    
    create_list(list, ranks, disps, sizes, &ranks_size, &total_size,
    temp_list, size, comm_size, my_rank);
    
    set_halo_list halo_list= (set_halo_list) malloc(sizeof(set_halo_list_core));
    halo_list->set = set;
    halo_list->size = total_size;
    halo_list->ranks = ranks;
    halo_list->ranks_size = ranks_size;
    halo_list->disps = disps;
    halo_list->sizes = sizes;
    halo_list->list = list;
    
    OP_export_sets_list[set->index] = halo_list;
    
}

void create_nonexec_set_import_list(op_set set, int* temp_list, int size, int comm_size, int my_rank)
{
    int* ranks = (int *) malloc(comm_size*sizeof(int));
    int* list = (int *) malloc((size/2)*sizeof(int));
    int* disps = (int *) malloc(comm_size*sizeof(int));
    int* sizes = (int *) malloc(comm_size*sizeof(int));
    
    int ranks_size = 0;
    int total_size = 0;  
    
    create_list(list, ranks, disps, sizes, &ranks_size, &total_size,
    temp_list, size, comm_size, my_rank);
    
    set_halo_list halo_list= (set_halo_list) malloc(sizeof(set_halo_list_core));
    halo_list->set = set;
    halo_list->size = total_size;
    halo_list->ranks = ranks;
    halo_list->ranks_size = ranks_size;
    halo_list->disps = disps;
    halo_list->sizes = sizes;
    halo_list->list = list;
    
    OP_import_nonexec_sets_list[set->index] = halo_list;
    
}

void create_map_export_list(op_map map, int* temp_list, int size, int comm_size, int my_rank)
{
    int* ranks = (int *) malloc(comm_size*sizeof(int));
    int* list = (int *) malloc((size/2)*sizeof(int));
    int* disps = (int *) malloc(comm_size*sizeof(int));
    int* sizes = (int *) malloc(comm_size*sizeof(int));
    
    int ranks_size = 0;
    int total_size = 0;  
    
    create_list(list, ranks, disps, sizes, &ranks_size, &total_size,
    temp_list, size, comm_size, my_rank);
    
    map_halo_list halo_list= (map_halo_list) malloc(sizeof(map_halo_list_core));
    halo_list->map = map;
    halo_list->size = total_size;
    halo_list->ranks = ranks;
    halo_list->ranks_size = ranks_size;
    halo_list->disps = disps;
    halo_list->sizes = sizes;
    halo_list->list = list;
    
    OP_export_maps_list[map->index] = halo_list;
    
}




void create_map_import_list(op_map map, int* temp_list, int total_size, int* ranks, 
    int* sizes, int ranks_size, int comm_size, int my_rank)
{
    int* disps = (int *) malloc(comm_size*sizeof(int));
    disps[0] = 0;
    for(int i=0; i<ranks_size; i++)
    {
    	if(i>0)disps[i] = disps[i-1]+sizes[i-1]; 	
    }
    
    map_halo_list halo_list= (map_halo_list) malloc(sizeof(map_halo_list_core));
    halo_list->map = map;
    halo_list->size = total_size;
    halo_list->ranks = ranks;
    halo_list->ranks_size = ranks_size;
    halo_list->disps = disps;
    halo_list->sizes = sizes;
    halo_list->list = temp_list;
    
    OP_import_maps_list[map->index] = halo_list;
    
}

void create_set_import_list(op_set set, int* temp_list, int total_size, int* ranks, 
    int* sizes, int ranks_size, int comm_size, int my_rank)
{
    int* disps = (int *) malloc(comm_size*sizeof(int));
    disps[0] = 0;
    for(int i=0; i<ranks_size; i++)
    {
    	if(i>0)disps[i] = disps[i-1]+sizes[i-1]; 	
    }
    
    set_halo_list halo_list= (set_halo_list) malloc(sizeof(set_halo_list_core));
    halo_list->set = set;
    halo_list->size = total_size;
    halo_list->ranks = ranks;
    halo_list->ranks_size = ranks_size;
    halo_list->disps = disps;
    halo_list->sizes = sizes;
    halo_list->list = temp_list;
    
    OP_import_sets_list[set->index] = halo_list;
    
}

void create_nonexec_set_export_list(op_set set, int* temp_list, int total_size, int* ranks, 
    int* sizes, int ranks_size, int comm_size, int my_rank)
{
    int* disps = (int *) malloc(comm_size*sizeof(int));
    disps[0] = 0;
    for(int i=0; i<ranks_size; i++)
    {
    	if(i>0)disps[i] = disps[i-1]+sizes[i-1]; 	
    }
    
    set_halo_list halo_list= (set_halo_list) malloc(sizeof(set_halo_list_core));
    halo_list->set = set;
    halo_list->size = total_size;
    halo_list->ranks = ranks;
    halo_list->ranks_size = ranks_size;
    halo_list->disps = disps;
    halo_list->sizes = sizes;
    halo_list->list = temp_list;
    
    OP_export_nonexec_sets_list[set->index] = halo_list;
}


/**--------------------------- Halo List Creation ---------------------------**/
void op_halo_create()
{
    //declare timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    double time;
    double max_time;
    op_timers(&cpu_t1, &wall_t1); //timer start for list create
    
    //create new communicator for OP mpi operation
    int my_rank, comm_size;
    MPI_Comm_dup(MPI_COMM_WORLD, &OP_MPI_WORLD);
    MPI_Comm_rank(OP_MPI_WORLD, &my_rank);
    MPI_Comm_size(OP_MPI_WORLD, &comm_size);
    
    /* Compute global partition range information for each set*/
    int** part_range = (int **) malloc(OP_set_index*sizeof(int*));
    if(part_range == NULL) {
    	printf(" op_list_create -- error allocating memory: int** part_range\n");
    	exit(-1);
    }
    get_part_range(part_range,my_rank,comm_size, OP_MPI_WORLD);
    
    OP_export_sets_list = (set_halo_list *)malloc(OP_set_index*sizeof(set_halo_list));
    OP_export_maps_list = (map_halo_list *)malloc(OP_map_index*sizeof(map_halo_list));
    
    
    /*----- STEP 1 - construct export lists for mappings and execute sets-----*/
    
    //declare temporaty scratch variables to hold set export lists and mapping 
    //table export lists
    int s_i;
    int* set_list;
  
    int* map_list[OP_map_index];
    int cap_s = 1000; //keep track of the temp array capacities
    
    
    for(int s=0; s<OP_set_index; s++){ //for each set
    	op_set set=OP_set_list[s];
    	
    	//create a temporaty scratch space to hold export list for this set
    	s_i = 0;cap_s = 1000;
    	set_list = (int *)malloc(cap_s*sizeof(int));
    	if(set_list == NULL) {
    	    printf(" op_list_create -- error allocating memory: int* set_list\n");
    	    exit(-1);
    	}
    	
    	for(int e=0; e<set->size;e++){//for each elment of this set
    	    for(int m=0; m<OP_map_index; m++) { //for each maping table
    	    	op_map map=OP_map_list[m];
    	    	
    	    	if(compare_sets(map->from,set)==1) //need to select mappings FROM this set
    	    	{
    	    	    int part, local_index;
    	    	    for(int j=0; j<map->dim; j++) { //for each element pointed at by this entry
    	    	    	part = get_partition(map->map[e*map->dim+j],
    	    	    	    part_range[map->to->index],&local_index,comm_size);
    	    	    	if(s_i>=cap_s)
    	    	    	{
    	    	    	    cap_s = cap_s*2;
    	    	    	    set_list = (int *)realloc(set_list,cap_s*sizeof(int));
    	    	    	    if(set_list == NULL) {
    	    	    	    	printf(" op_list_create -- error reallocating memory: int* set_list\n");
    	    	    	    	exit(-1);
    	    	    	    }
    	    	    	}
    	    	    	
    	    	    	if(part != my_rank){
    	    	    	    set_list[s_i++] = part; //add to set export list
    	    	    	    set_list[s_i++] = e;
    	    	    	}
    	    	    }
    	    	}
    	    }
    	}
    	
    	//create set export list
    	//printf("creating set export list for set %10s of size %d\n",set->name,s_i);
    	create_set_export_list(set,set_list,s_i, comm_size, my_rank);
    	for(int m=0; m<OP_map_index; m++) { //for each maping table create export list
    	    op_map map=OP_map_list[m];
    	    if(compare_sets(map->from,set)==1) //need to select mappings FROM this set
    	    {
    	    	//create map export list: the union of all entries in mapping tables that is FROM this set
    	    	//printf("creating map export list for map %10s of size %d\n",map->name,s_i);
    	    	create_map_export_list(map,set_list, s_i, comm_size, my_rank);    	    	
    	    }
    	}
    	free(set_list);//free temp list
    	
    }
    
    
    /*for(int s=0; s<OP_set_index; s++){ //for each set
    	op_set set=OP_set_list[s];
    	print_set_list(set, my_rank, 1, 0);
    }*/
    
    
    /*---- STEP 2 - construct import lists for mappings and execute sets------*/
    
    OP_import_sets_list = (set_halo_list *)malloc(OP_set_index*sizeof(set_halo_list));
    OP_import_maps_list = (map_halo_list *)malloc(OP_map_index*sizeof(map_halo_list));
    
    /* Mappings first................*/
    int *neighbors, *sizes;
    int ranks_size;
    
    for(int m=0; m<OP_map_index; m++) { //for each maping table
    	op_map map=OP_map_list[m];
    	
    	//-----Discover neighbors-----
    	ranks_size = 0;
    	neighbors = malloc(comm_size*sizeof(int));
    	if(neighbors == NULL) {
    	    printf(" op_list_create -- error allocating memory: int *neighbors\n");
    	    exit(-1);
    	}
    	sizes = malloc(comm_size*sizeof(int));
    	if(sizes == NULL) {
    	    printf(" op_list_create -- error allocating memory: int *sizes\n");
    	    exit(-1);
    	}
    	
    	map_halo_list list=OP_export_maps_list[map->index];
    	
    	find_neighbors_map(list, neighbors, sizes, &ranks_size, my_rank, comm_size, OP_MPI_WORLD);
    	MPI_Request request_send[list->ranks_size];
    	
    	int* rbuf, cap = 0, index = 0;
    	
    	for(int i=0; i< ranks_size; i++) cap = cap + sizes[i];    	     
    	int* temp = malloc(cap*sizeof(int));
    	if(temp == NULL) {
    	    printf(" op_list_create -- error allocating memory: int* temp\n");
    	    exit(-1);
    	}
 
    	for(int i=0; i<list->ranks_size; i++) {
    	    //printf("export map %10s to %d from rank %d, list of size %d \n",
    	    //    map->name,list->ranks[i],my_rank,list->sizes[i]);
    	    int* sbuf = &list->list[list->disps[i]];
    	    MPI_Isend( sbuf,  list->sizes[i],  MPI_INT, list->ranks[i], m,
    	    	OP_MPI_WORLD, &request_send[i]);
    	}
     
    	//import this list from those neighbors
    	for(int i=0; i<ranks_size; i++) {
    	    //printf("on rank %d, import map %10s from %d list of size %d \n",
    	    //    my_rank,map->name,neighbors[i],sizes[i]);
    	    rbuf = malloc(sizes[i]*sizeof(int));
    	    if(rbuf == NULL) {
    	    	printf(" op_list_create -- error allocating memory: int* rbuf\n");
    	    	exit(-1);
    	    }
    	
    	    MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], m, OP_MPI_WORLD,
    	    	MPI_STATUSES_IGNORE );
    	    memcpy(&temp[index],(void *)&rbuf[0],sizes[i]*sizeof(int));
    	    index = index + sizes[i];
    	    free(rbuf);
    	}
 
    	MPI_Waitall(list->ranks_size,request_send, MPI_STATUSES_IGNORE );
    	
    	//create import lists
    	//printf("creating importlist of with number of neighbors %d\n",ranks_size);
    	create_map_import_list(map, temp, index,neighbors, sizes, 
    	    ranks_size, comm_size, my_rank);
    	
    }
    
        
    /* Sets next..............*/
    
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=OP_set_list[s];
    	
    	//-----Discover neighbors-----
    	ranks_size = 0;
    	neighbors = malloc(comm_size*sizeof(int));
    	if(neighbors == NULL) {
    	    printf(" op_list_create -- error allocating memory: int *neighbors\n");
    	    exit(-1);
    	}
    	sizes = malloc(comm_size*sizeof(int));
    	if(sizes == NULL) {
    	    printf(" op_list_create -- error allocating memory: int *sizes\n");
    	    exit(-1);
    	}
    	
    	set_halo_list list = OP_export_sets_list[set->index];
    	
    	find_neighbors_set(list,neighbors,sizes,&ranks_size,my_rank,comm_size, OP_MPI_WORLD);
    	MPI_Request request_send[list->ranks_size];
    	
    	int* rbuf, cap = 0, index = 0;
    	
    	for(int i=0; i<list->ranks_size; i++) {
    	    //printf("export from %d to %d set %10s, list of size %d \n",
    	    //my_rank,list->ranks[i],set->name,list->sizes[i]);
    	    int* sbuf = &list->list[list->disps[i]];
    	    MPI_Isend( sbuf,  list->sizes[i],  MPI_INT, list->ranks[i], s,
    	    	OP_MPI_WORLD, &request_send[i] );
    	}
    	
    	for(int i=0; i< ranks_size; i++) cap = cap + sizes[i];
    	int* temp = malloc(cap*sizeof(int));
    	if(temp == NULL) {
    	    printf(" op_list_create -- error allocating memory: int* temp\n");
    	    exit(-1);
    	}
    	
    	//import this list from those neighbors
    	for(int i=0; i<ranks_size; i++) {
    	    //printf("import from %d to %d set %10s, list of size %d\n",
    	    //neighbors[i], my_rank, set->name, sizes[i]);
    	    rbuf = malloc(sizes[i]*sizeof(int));
    	    if(rbuf == NULL) {
    	    	printf(" op_list_create -- error allocating memory: int* rbuf\n");
    	    	exit(-1);
    	    }
    	    MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i],s, OP_MPI_WORLD,
    	    	MPI_STATUSES_IGNORE );
    	    memcpy(&temp[index],(void *)&rbuf[0],sizes[i]*sizeof(int));
    	    index = index + sizes[i];
    	    free(rbuf);
    	}
    	
    	MPI_Waitall(list->ranks_size,request_send, MPI_STATUSES_IGNORE );
    	
    	//create import lists
    	//printf("creating importlist with number of neighbors %d\n",ranks_size);
    	create_set_import_list(set, temp, index,neighbors, sizes, 
    	    ranks_size, comm_size, my_rank);
    }
    
    
    /*-- STEP 3 - Exchange mapping table entries using the import/export lists--*/
  
    for(int m=0; m<OP_map_index; m++) { //for each maping table
    	op_map map=OP_map_list[m];
    	map_halo_list i_list = OP_import_maps_list[map->index];
    	map_halo_list e_list = OP_export_maps_list[map->index];
  	  
    	MPI_Request request_send[e_list->ranks_size];

    	//prepare bits of the mapping tables to be exported
    	int** sbuf = (int **)malloc(e_list->ranks_size*sizeof(int *));
    	if(sbuf == NULL) {
    	    	printf(" op_list_create -- error allocating memory: int** sbuf\n");
    	    	exit(-1);
    	}
    	
    	for(int i=0; i < e_list->ranks_size; i++) {
    	    sbuf[i] = (int *)malloc(e_list->sizes[i]*map->dim*sizeof(int));
    	    if(sbuf[i] == NULL) {
    	    	printf(" op_list_create -- error allocating memory: int* sbuf[i]\n");
    	    	exit(-1);
    	    }
    	    
    	    for(int j = 0; j < e_list->sizes[i]; j++)
    	    {
    	    	for(int p = 0; p < map->dim; p++)
    	    	{
    	    	    sbuf[i][j*map->dim+p] =
    	    	    map->map[map->dim*(e_list->list[e_list->disps[i]+j])+p];
    	    	}
    	    }
    	    //printf("\n export from %d to %d map %10s, number of elements of size %d | sending:\n ",
    	    //    my_rank,e_list.ranks[i],map.name,e_list.sizes[i]);
    	    MPI_Isend(sbuf[i],  map->dim*e_list->sizes[i],  MPI_INT, e_list->ranks[i],
    	    	m, OP_MPI_WORLD, &request_send[i]);
    	}
      
    	//prepare space for the incomming mapping tables - realloc each
    	//mapping tables in each mpi process
    	OP_map_list[map->index]->map = realloc(OP_map_list[map->index]->map,
    	    (map->dim*(map->from->size+i_list->size))*sizeof(int));
        if(OP_map_list[map->index]->map == NULL) {
    	    	printf(" op_list_create -- error reallocating memory: OP_map_list[map->index]->map\n");
    	    	exit(-1);
    	}
    	    
    	int init = map->dim*(map->from->size);
    	for(int i=0; i<i_list->ranks_size; i++) {
    	    //printf("\n imported on to %d map %10s, number of elements of size %d | recieving: ",
    	    //	  my_rank, map->name, i_list->size);
    	    MPI_Recv(&(OP_map_list[map->index]->map[init+i_list->disps[i]*map->dim]),
    	    	map->dim*i_list->sizes[i], MPI_INT, i_list->ranks[i], m,
    	    	OP_MPI_WORLD, MPI_STATUSES_IGNORE);
    	}

    	MPI_Waitall(e_list->ranks_size,request_send, MPI_STATUSES_IGNORE );
    	for(int i=0; i < e_list->ranks_size; i++) free(sbuf[i]); free(sbuf);    	
    }
    
    
    
    /*-- STEP 4 - Create import lists for non-execute set elements using mapping
    table entries including the additional mapping table entries --*/
    
    OP_import_nonexec_sets_list = (set_halo_list *)malloc(OP_set_index*sizeof(set_halo_list));
    OP_export_nonexec_sets_list = (set_halo_list *)malloc(OP_set_index*sizeof(set_halo_list));
     
    //declare temporaty scratch variables to hold non-exec set export lists
    s_i = 0;
    set_list = NULL;  
    cap_s = 1000; //keep track of the temp array capacity
  
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=OP_set_list[s];
    	set_halo_list exec_set_list=OP_import_sets_list[set->index];   
      
    	//create a temporaty scratch space to hold nonexec export list for this set
    	s_i = 0;
    	set_list = malloc(cap_s*sizeof(int));     
        if(set_list == NULL) {
            printf(" op_list_create -- error allocating memory: set_list\n");
            exit(-1);
        }
    	
    	for(int m=0; m<OP_map_index; m++) { //for each maping table
    	    op_map map=OP_map_list[m];
    	    map_halo_list exec_map_list=OP_import_maps_list[map->index];
 
    	    if(compare_sets(map->to,set)==1) //need to select mappings TO this set
    	    {
    	    	//for each entry in this mapping table: original+execlist
    	    	int len = map->from->size+exec_map_list->size;
    	    	for(int e = 0; e<len; e++)
    	    	{
    	    	    int part;
    	    	    int local_index;
    	    	    for(int j=0; j < map->dim; j++) { //for each element pointed at by this entry
    	    	    	part = get_partition(map->map[e*map->dim+j],
    	    	    	    part_range[map->to->index],&local_index,comm_size);
              	      
    	    	    	if(s_i>=cap_s)
    	    	    	{
    	    	    	    cap_s = cap_s*2;
    	    	    	    set_list = realloc(set_list,cap_s*sizeof(int));
    	    	    	    if(set_list == NULL) {
    	    	    	    	printf(" op_list_create -- error reallocating memory: set_list\n");
    	    	    	    	exit(-1);
    	    	    	    }
    	    	    	}
              	      
    	    	    	if(part != my_rank)
    	    	    	{
    	    	    	    int found = -1;
    	    	    	    //check in exec list
    	    	    	    int rank = binary_search(exec_set_list->ranks,
    	    	    	    	part, 0, exec_set_list->ranks_size-1);
              	      	  
    	    	    	    if(rank >= 0)
    	    	    	    {
    	    	    	    	found = binary_search(exec_set_list->list,
    	    	    	    	    local_index, exec_set_list->disps[rank],
    	    	    	    	    exec_set_list->disps[rank]+exec_set_list->sizes[rank]-1);
    	    	    	    }
              	      	  
    	    	    	    if(found < 0){
    	    	    	    	// not in this partition and not found in exec list
    	    	    	    	//add to non-execute set_list
    	    	    	    	set_list[s_i++] = part;
    	    	    	    	set_list[s_i++] = local_index;
    	    	    	    }
    	    	    	}
    	    	    }
    	    	}
    	    }
    	}
      
    	//create non-exec set import list
    	//printf("creating non-exec import list of size %d\n",s_i);
    	create_nonexec_set_import_list(set,set_list, s_i, comm_size, my_rank);
    	free(set_list);//free temp list
      
    }
    
    
    /*----------- STEP 5 - construct non-execute set export lists -------------*/
  
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=OP_set_list[s];
      
    	//-----Discover neighbors-----
    	ranks_size = 0;
    	neighbors = malloc(comm_size*sizeof(int));
    	if(neighbors == NULL) {
    	    printf(" op_list_create -- error allocating memory: int *neighbors\n");
    	    exit(-1);
    	}
    	sizes = malloc(comm_size*sizeof(int));
    	if(sizes == NULL) {
    	    printf(" op_list_create -- error allocating memory: int *sizes\n");
    	    exit(-1);
    	}
      
    	set_halo_list list=OP_import_nonexec_sets_list[set->index];
    	find_neighbors_set(list,neighbors,sizes,&ranks_size,my_rank,comm_size, OP_MPI_WORLD);
            
    	MPI_Request request_send[list->ranks_size];
            
    	int* rbuf, cap = 0, index = 0;
      
    	for(int i=0; i<list->ranks_size; i++) {
    	    //printf("import to %d from %d set %10s, nonexec list of size %d | sending:\n",
    	    //    my_rank,list->ranks[i],set->name,list->sizes[i]);
    	    int* sbuf = &list->list[list->disps[i]];
    	    MPI_Isend( sbuf,  list->sizes[i],  MPI_INT, list->ranks[i], s,
    	    	OP_MPI_WORLD, &request_send[i] );  
    	}
      
    	for(int i=0; i< ranks_size; i++) cap = cap + sizes[i];     	     
    	int* temp = malloc(cap*sizeof(int));
    	if(temp == NULL) {
    	    printf(" op_list_create -- error allocating memory: int* temp\n");
    	    exit(-1);
    	}
     
    	//export this list to those neighbors
    	for(int i=0; i<ranks_size; i++) {
    	    //printf("export to %d from %d set %10s, list of size %d | recieving:\n",
    	    //    neighbors[i], my_rank, set->name, sizes[i]);
    	    rbuf = malloc(sizes[i]*sizeof(int));
    	    if(sizes == NULL) {
    	    	printf(" op_list_create -- error allocating memory: rbuf\n");
    	    	exit(-1);
    	    }
    	    MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i],s, OP_MPI_WORLD,
    	    	MPI_STATUSES_IGNORE );
    	    memcpy(&temp[index],(void *)&rbuf[0],sizes[i]*sizeof(int));
    	    index = index + sizes[i];
    	    free(rbuf);
    	}
      
    	MPI_Waitall(list->ranks_size,request_send, MPI_STATUSES_IGNORE );
     
    	//create import lists
    	//printf("creating nonexec set export list with number of neighbors %d\n",ranks_size);
    	create_nonexec_set_export_list(set, temp, index, neighbors, sizes, 
    	    ranks_size, comm_size, my_rank);
    }
    
    
    /*-STEP 6 - Exchange execute set elements/data using the import/export lists--*/
    
    for(int s=0; s<OP_set_index; s++){ //for each set
    	op_set set=OP_set_list[s];
    	set_halo_list i_list = OP_import_sets_list[set->index];
    	set_halo_list e_list = OP_export_sets_list[set->index];
      
    	//for each data array
    	for(int d=0; d<OP_dat_index; d++){
    	    op_dat dat=OP_dat_list[d];
      	  
    	    if(compare_sets(set,dat->set)==1)//if this data array is defined on this set
    	    {
    	    	//printf("on rank %d, The data array is %10s\n",my_rank,dat->name);
    	    	MPI_Request request_send[e_list->ranks_size];
      	      
    	    	//prepare execute set element data to be exported
    	    	char** sbuf = malloc(e_list->ranks_size*sizeof(char *));
    	    	if(sbuf == NULL) {
    	    	    printf(" op_list_create -- error allocating memory: char** sbuf\n");
    	    	    exit(-1);
    	    	}
    	    
    	    	for(int i=0; i < e_list->ranks_size; i++) {
    	    	    sbuf[i] = malloc(e_list->sizes[i]*dat->size);
    	    	    if(sbuf == NULL) {
    	    	    	printf(" op_list_create -- error allocating memory: sbuf[i]\n");
    	    	    	exit(-1);
    	    	    }
    	    	    for(int j = 0; j < e_list->sizes[i]; j++)
    	    	    {
    	    	    	int set_elem_index = e_list->list[e_list->disps[i]+j];
    	    	    	memcpy(&sbuf[i][j*dat->size],(void *)&dat->data[dat->size*(set_elem_index)],dat->size);
    	    	    }
    	    	    //printf("export from %d to %d data %10s, number of elements of size %d | sending:\n ",
    	    	    //    my_rank,e_list->ranks[i],dat->name,e_list->sizes[i]);
    	    	    MPI_Isend(sbuf[i],  dat->size*e_list->sizes[i],  MPI_CHAR, e_list->ranks[i],
    	    	    	d, OP_MPI_WORLD, &request_send[i]);
    	    	}
      	      
    	    	//prepare space for the incomming data - realloc each
    	    	//data array in each mpi process
    	    	OP_dat_list[dat->index]->data = realloc(OP_dat_list[dat->index]->data,
    	    	    (set->size+i_list->size)*dat->size);
    	    	if(OP_dat_list[dat->index]->data == NULL) {
    	    	    printf(" op_list_create -- error reallocating memory: \
    	    	    	OP_dat_list[dat->index]->data\n");
    	    	    exit(-1);
    	    	}
    	    	
    	    	int init = set->size*dat->size;
    	    	for(int i=0; i<i_list->ranks_size; i++) {
    	    	    MPI_Recv(&(OP_dat_list[dat->index]->data[init+i_list->disps[i]*dat->size]),
    	    	    	dat->size*i_list->sizes[i], MPI_CHAR, i_list->ranks[i], d,
    	    	    	OP_MPI_WORLD, MPI_STATUSES_IGNORE);
    	    	}
      	      
    	    	MPI_Waitall(e_list->ranks_size,request_send, MPI_STATUSES_IGNORE );
    	    	for(int i=0; i<e_list->ranks_size; i++) free(sbuf[i]); free(sbuf);
    	    	//printf("imported on to %d data %10s, number of elements of size %d | recieving:\n ",
    	    	//	  my_rank, dat->name, i_list->size);
    	    }
    	}
    }
    
    /*-STEP 7 - Exchange non-execute set elements/data using the import/export lists--*/

    for(int s=0; s<OP_set_index; s++){ //for each set
    	op_set set=OP_set_list[s];
    	set_halo_list i_list = OP_import_nonexec_sets_list[set->index];
    	set_halo_list e_list = OP_export_nonexec_sets_list[set->index];
        
    	//for each data array
    	for(int d=0; d<OP_dat_index; d++){
    	    op_dat dat=OP_dat_list[d];
      	  
    	    if(compare_sets(set,dat->set)==1)//if this data array is defined on this set
    	    {
    	    	//printf("on rank %d, The data array is %10s\n",my_rank,dat->name);
    	    	MPI_Request request_send[e_list->ranks_size];
      	      
    	    	//prepare non-execute set element data to be exported
    	    	char** sbuf = malloc(e_list->ranks_size*sizeof(char *));
    	    	if(sbuf == NULL) {
    	    	    printf(" op_list_create -- error allocating memory: char** sbuf\n");
    	    	    exit(-1);
    	    	}
    	    	
    	    	for(int i=0; i < e_list->ranks_size; i++) {
    	    	    sbuf[i] = malloc(e_list->sizes[i]*dat->size);
    	    	    if(sbuf[i] == NULL) {
    	    	    	printf(" op_list_create -- error allocating memory: char* sbuf[i]\n");
    	    	    	exit(-1);
    	    	    }
    	    	
    	    	    for(int j = 0; j < e_list->sizes[i]; j++)
    	    	    {
    	    	    	int set_elem_index = e_list->list[e_list->disps[i]+j];
    	    	    	memcpy(&sbuf[i][j*dat->size],
    	    	    	    (void *)&dat->data[dat->size*(set_elem_index)],dat->size);
    	    	    }
    	    	    MPI_Isend(sbuf[i],  dat->size*e_list->sizes[i],  MPI_CHAR, e_list->ranks[i],
    	    	    	d, OP_MPI_WORLD, &request_send[i]);
    	    	}
      	      
    	    	//prepare space for the incomming nonexec-data - realloc each
    	    	//data array in each mpi process
    	    	set_halo_list exec_i_list = OP_import_sets_list[set->index];      
      
    	    	OP_dat_list[dat->index]->data = realloc(OP_dat_list[dat->index]->data,
    	    	    (set->size+exec_i_list->size+i_list->size)*dat->size);
    	    	if(OP_dat_list[dat->index]->data == NULL) {
    	    	    printf(" op_list_create -- error reallocating memory: \
    	    	    	OP_dat_list[dat->index]->data\n");
    	    	    exit(-1);
    	    	}
    	    	    
    	    	int init = (set->size+exec_i_list->size)*dat->size;
      	      
    	    	for(int i=0; i < i_list->ranks_size; i++) {
    	    	    MPI_Recv(&(OP_dat_list[dat->index]->data[init+i_list->disps[i]*dat->size]),
    	    	    	dat->size*i_list->sizes[i], MPI_CHAR, i_list->ranks[i], d,
    	    	    	OP_MPI_WORLD, MPI_STATUSES_IGNORE);
    	    	}
      	      
    	    	MPI_Waitall(e_list->ranks_size,request_send, MPI_STATUSES_IGNORE );
    	    	for(int i=0; i < e_list->ranks_size; i++) free(sbuf[i]); free(sbuf);
    	    }
    	}
    }
    
    
    
    /*-STEP 8 ----------------- Renumber Mapping tables-----------------------*/
	
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=OP_set_list[s];
      
    	for(int m=0; m<OP_map_index; m++) { //for each maping table
    	    op_map map=OP_map_list[m];
      	        	 
    	    if(compare_sets(map->to,set)==1) //need to select mappings TO this set
    	    {	
    	    	set_halo_list exec_set_list=OP_import_sets_list[set->index];
    	    	set_halo_list nonexec_set_list=OP_import_nonexec_sets_list[set->index];
      	  
    	    	map_halo_list exec_map_list=OP_import_maps_list[map->index];
              
    	    	//for each entry in this mapping table: original+execlist
    	    	int len = map->from->size+exec_map_list->size;
    	    	for(int e = 0; e < len; e++)
    	    	{
    	    	    for(int j=0; j < map->dim; j++) { //for each element pointed at by this entry
    	    	    	
    	    	    	int part;
    	    	    	int local_index = 0;
    	    	    	part = get_partition(map->map[e*map->dim+j],
    	    	    	    part_range[map->to->index],&local_index,comm_size);
              	                    	      
    	    	    	if(part == my_rank)
    	    	    	{
    	    	    	    OP_map_list[map->index]->map[e*map->dim+j] = local_index;
    	    	    	}
    	    	    	else
    	    	    	{
    	    	    	    int found = -1;
    	    	    	    //check in exec list
    	    	    	    int rank1 = binary_search(exec_set_list->ranks,
    	    	    	    	part, 0, exec_set_list->ranks_size-1);
    	    	    	    //check in nonexec list
    	    	    	    int rank2 = binary_search(nonexec_set_list->ranks,
    	    	    	    	part, 0, nonexec_set_list->ranks_size-1);
              	      	  
    	    	    	    if(rank1 >=0)
    	    	    	    {
    	    	    	    	found = binary_search(exec_set_list->list,
    	    	    	    	    local_index, exec_set_list->disps[rank1],
    	    	    	    	    exec_set_list->disps[rank1]+exec_set_list->sizes[rank1]-1);

              	      	      	if(found>=0) 
              	      	      	{
              	      	      	    OP_map_list[map->index]->map[e*map->dim+j] = 
              	      	      	    found + map->to->size ;
              	      	      	}
              	      	    }
              	      	  
              	      	    if(rank2 >=0 & found <0)
              	      	    {
              	      	    	found = binary_search(nonexec_set_list->list,
              	      	    	    local_index, nonexec_set_list->disps[rank2],
              	      	    	    nonexec_set_list->disps[rank2]+
              	      	    	    nonexec_set_list->sizes[rank2]-1);
              	      	    	
              	      	    	if(found>=0)
              	      	      	{
              	      	      	    OP_map_list[map->index]->map[e*map->dim+j] = 
              	      	      	    found + set->size + exec_set_list->size;
              	      	      	}
              	      	    }
              	      	  
              	      	    if(found < 0) 
              	      	    printf("ERROR: Set %10s Element %d needed on rank %d \
              	      	    	from partition %d\n",set->name, local_index, my_rank,part );
              	      	}
              	    }
              	}
            }
        }
    }
    
    
    //set dirty bits of all data arrays to 0
    //for each data array
    
    dirtybit = (int *)malloc(OP_dat_index*sizeof(int));
    if(dirtybit == NULL) {
    	printf(" op_list_create -- error allocating memory:  dirtybit\n"); 
    	exit(-1);
    }
    
    for(int d=0; d<OP_dat_index; d++){
    	op_dat dat=OP_dat_list[d];
    	
    	dirtybit[dat->index] = 0;
    }
    
    op_timers(&cpu_t2, &wall_t2);  //timer stop for list create    
    //compute import/export lists creation time
    time = wall_t2-wall_t1;
    MPI_Reduce(&time,&max_time,1,MPI_DOUBLE, MPI_MAX,0, OP_MPI_WORLD);
    
    
    //compute average halo size in Bytes
    int tot_halo_size = 0;
    for(int s = 0; s< OP_set_index; s++){
    	op_set set=OP_set_list[s];   
    	
    	for(int d=0; d<OP_dat_index; d++){
    	    op_dat dat=OP_dat_list[d];
    	    
    	     if(compare_sets(dat->set,set)==1)
    	     {
    	     	 set_halo_list exec_imp = OP_import_sets_list[set->index];
    	     	 set_halo_list nonexec_imp= OP_import_nonexec_sets_list[set->index];
    	     	 
    	     	 tot_halo_size = tot_halo_size + exec_imp->size + nonexec_imp->size;
    	     }
    	}
    }
    tot_halo_size = tot_halo_size * sizeof(double);
    int avg_halo_size;
    MPI_Reduce(&tot_halo_size,&avg_halo_size,1,MPI_INT, MPI_SUM,0, OP_MPI_WORLD);
    
    //print performance results
    if(my_rank==0)
    {
    	printf("Max total halo creation time = %lf\n",max_time);
    	printf("Average Halo size = %d Bytes\n",tot_halo_size/comm_size);
    }   
}



/**--------------------------- Clean-up Halo Lists --------------------------**/

void op_halo_destroy()
{
    for(int s = 0; s< OP_set_index; s++){
    	op_set set=OP_set_list[s]; 
    	    	
    	free(OP_import_sets_list[set->index]->ranks);
    	free(OP_import_sets_list[set->index]->disps);
    	free(OP_import_sets_list[set->index]->sizes);
    	free(OP_import_sets_list[set->index]->list);
    	
    	free(OP_import_nonexec_sets_list[set->index]->ranks);
    	free(OP_import_nonexec_sets_list[set->index]->disps);
    	free(OP_import_nonexec_sets_list[set->index]->sizes);
    	free(OP_import_nonexec_sets_list[set->index]->list);
    	
    	free(OP_export_sets_list[set->index]->ranks);
    	free(OP_export_sets_list[set->index]->disps);
    	free(OP_export_sets_list[set->index]->sizes);
    	free(OP_export_sets_list[set->index]->list);
    	
    	free(OP_export_nonexec_sets_list[set->index]->ranks);
    	free(OP_export_nonexec_sets_list[set->index]->disps);
    	free(OP_export_nonexec_sets_list[set->index]->sizes);
    	free(OP_export_nonexec_sets_list[set->index]->list);    	
    	
    }
    
    for(int m = 0; m< OP_map_index; m++){
    	op_map map=OP_map_list[m]; 
    	    	
    	free(OP_import_maps_list[map->index]->ranks);
    	free(OP_import_maps_list[map->index]->disps);
    	free(OP_import_maps_list[map->index]->sizes);
    	free(OP_import_maps_list[map->index]->list);
    	
    	free(OP_export_maps_list[map->index]->ranks);
    	free(OP_export_maps_list[map->index]->disps);
    	free(OP_export_maps_list[map->index]->sizes);
    	free(OP_export_maps_list[map->index]->list);    	
    }
    
    
}





/**-------------------------MPI Halo Exchange Functions----------------------**/

void exchange_halo(op_set set, op_arg arg)
{
    //if(arg.argtype == OP_ARG_DAT)
    //{
	op_dat dat = arg.dat;
	
	if((arg.idx != -1) && (arg.acc == OP_READ || arg.acc == OP_RW ) && 
	    (dirtybit[dat->index] == 1))
	{
	    //printf("Exchanging Halo of data array %10s\n",dat->name);
	    
	    set_halo_list imp_exec_list = OP_import_sets_list[dat->set->index];
	    set_halo_list imp_nonexec_list = OP_import_nonexec_sets_list[dat->set->index];
	    
	    set_halo_list exp_exec_list = OP_export_sets_list[dat->set->index];
	    set_halo_list exp_nonexec_list = OP_export_nonexec_sets_list[dat->set->index];
	    
	    //-------first exchange exec elements related to this data array--------
	
	    //sanity checks
	    if(compare_sets(imp_exec_list->set,dat->set)==0) 
		{ printf("Error: Import list and set mismatch\n"); exit(2);}
	    if(compare_sets(exp_exec_list->set,dat->set)==0) 
		{printf("Error: Export list and set mismatch\n"); exit(2);}
		    
	    MPI_Request request_send[exp_exec_list->ranks_size];
	    char* sbuf[exp_exec_list->ranks_size];
	    //prepare execute set element data to be exported
	    for(int i=0; i<exp_exec_list->ranks_size; i++) {
		sbuf[i] = malloc(exp_exec_list->sizes[i]*dat->size);
		if(sbuf[i] == NULL) {
		    printf(" exchange_halo -- error allocating memory:  sbuf[i]\n");
		    exit(-1);
		}
	
		for(int j = 0; j < exp_exec_list->sizes[i]; j++)
		{
		    int set_elem_index = exp_exec_list->list[exp_exec_list->disps[i]+j];
		    memcpy(&sbuf[i][j*dat->size],
			(void *)&dat->data[dat->size*(set_elem_index)],dat->size);
		}
		//printf("export from %d to %d data %10s, number of elements of size %d | sending:\n ",
		//  	      my_rank, exp_exec_list->ranks[i], dat->name,exp_exec_list->sizes[i]);
		MPI_Isend(sbuf[i],  dat->size*exp_exec_list->sizes[i],  
		    MPI_CHAR, exp_exec_list->ranks[i],
		    dat->index, OP_MPI_WORLD, &request_send[i]);
	    }
	    
	    int init = dat->set->size*dat->size;
	    
	    for(int i=0; i < imp_exec_list->ranks_size; i++) {
		//printf("import on to %d from %d data %10s, number of elements of size %d | recieving:\n ",
		//  	  my_rank, imp_exec_list.ranks[i], dat.name, imp_exec_list.sizes[i]);
		MPI_Recv(&(OP_dat_list[dat->index]->data[init+imp_exec_list->disps[i]*dat->size]),
		    dat->size*imp_exec_list->sizes[i], MPI_CHAR, 
		    imp_exec_list->ranks[i], dat->index,
		    OP_MPI_WORLD, MPI_STATUSES_IGNORE);
	    }
	    MPI_Waitall(exp_exec_list->ranks_size,request_send, MPI_STATUSES_IGNORE );
	    for(int i=0; i<exp_exec_list->ranks_size; i++) free(sbuf[i]);
		
	
	    //-----second exchange nonexec elements related to this data array------
	    //sanity checks
	    if(compare_sets(imp_nonexec_list->set,dat->set)==0) 
		{ printf("Error: Non-Import list and set mismatch"); exit(2);}
	    if(compare_sets(exp_nonexec_list->set,dat->set)==0) 
		{printf("Error: Non-Export list and set mismatch"); exit(2);}
	
	    
	    MPI_Request request_nonexec_send[exp_nonexec_list->ranks_size];
	    char* sbuf_nonexec[exp_nonexec_list->ranks_size];
	    //prepare execute set element data to be exported
	    for(int i=0; i<exp_nonexec_list->ranks_size; i++) {
		sbuf_nonexec[i] = malloc(exp_nonexec_list->sizes[i]*dat->size);    	    
		if(sbuf_nonexec[i] == NULL) {
		    printf(" exchange_halo -- error allocating memory:  sbuf_nonexec[i]\n");
		    exit(-1);
		}
		
		for(int j = 0; j < exp_nonexec_list->sizes[i]; j++)
		{
		    int set_elem_index = exp_nonexec_list->list[exp_nonexec_list->disps[i]+j];
		    memcpy(&sbuf_nonexec[i][j*dat->size],
			(void *)&dat->data[dat->size*(set_elem_index)],dat->size);
		}
		
		MPI_Isend(sbuf_nonexec[i],  dat->size*exp_nonexec_list->sizes[i],  
		    MPI_CHAR, exp_nonexec_list->ranks[i],
		    dat->index, OP_MPI_WORLD, &request_nonexec_send[i]);
	    }
	    
	    int nonexec_init = (dat->set->size+imp_exec_list->size)*dat->size;
	    
	    for(int i=0; i<imp_nonexec_list->ranks_size; i++) {
		MPI_Recv(&(OP_dat_list[dat->index]->data[nonexec_init+imp_nonexec_list->disps[i]*dat->size]),
		    dat->size*imp_nonexec_list->sizes[i], MPI_CHAR, imp_nonexec_list->ranks[i], dat->index,
		    OP_MPI_WORLD, MPI_STATUSES_IGNORE);
	    }
	    MPI_Waitall(exp_nonexec_list->ranks_size,request_nonexec_send, MPI_STATUSES_IGNORE );
	    for(int i=0; i<exp_nonexec_list->ranks_size; i++) free(sbuf_nonexec[i]);
		    
	    //clear dirty bit
	    dirtybit[dat->index] = 0;
	}
   // }
    
}


void set_dirtybit(op_arg arg)
{
    //if(arg.argtype == OP_ARG_DAT)
    	if(arg.acc == OP_INC || arg.acc == OP_WRITE || arg.acc == OP_RW)
    	dirtybit[arg.dat->index] = 1;
}


void global_reduce(op_arg *arg)
{
    if(strcmp("double",arg->type)==0)
    {
	double result;
	if(arg->acc == OP_INC)//global reduction
	{
	    MPI_Reduce((double *)arg->data,&result,1,MPI_DOUBLE, MPI_SUM,0, OP_MPI_WORLD);
	    memcpy(arg->data, &result, sizeof(double));
	}
	else if(arg->acc == OP_MAX)//global maximum
	{
	    MPI_Reduce((double *)arg->data,&result,1,MPI_DOUBLE, MPI_MAX,0, OP_MPI_WORLD);
	    memcpy(arg->data, &result, sizeof(double));;
	}
	else if(arg->acc == OP_MIN)//global minimum              
	{
	    MPI_Reduce((double *)arg->data,&result,1,MPI_DOUBLE, MPI_MIN,0, OP_MPI_WORLD);
	    memcpy(arg->data, &result, sizeof(double));
	}
    }
    else if(strcmp("float",arg->type)==0)
    {
	float result;
	if(arg->acc == OP_INC)//global reduction
	{
	    MPI_Reduce((float *)arg->data,&result,1,MPI_FLOAT, MPI_SUM,0, OP_MPI_WORLD);
	    memcpy(arg->data, &result, sizeof(float));
	}
	else if(arg->acc == OP_MAX)//global maximum
	{
	    MPI_Reduce((float *)arg->data,&result,1,MPI_FLOAT, MPI_MAX,0, OP_MPI_WORLD);
	    memcpy(arg->data, &result, sizeof(float));;
	}
	else if(arg->acc == OP_MIN)//global minimum              
	{
	    MPI_Reduce((float *)arg->data,&result,1,MPI_FLOAT, MPI_MIN,0, OP_MPI_WORLD);
	    memcpy(arg->data, &result, sizeof(float));
	}
    }
    else if(strcmp("int",arg->type)==0)
    {
	int result;
	if(arg->acc == OP_INC)//global reduction
	{
	    MPI_Reduce((int *)arg->data,&result,1,MPI_INT, MPI_SUM,0, OP_MPI_WORLD);
	    memcpy(arg->data, &result, sizeof(int));
	}
	else if(arg->acc == OP_MAX)//global maximum
	{
	    MPI_Reduce((int *)arg->data,&result,1,MPI_INT, MPI_MAX,0, OP_MPI_WORLD);
	    memcpy(arg->data, &result, sizeof(int));;
	}
	else if(arg->acc == OP_MIN)//global minimum              
	{
	    MPI_Reduce((int *)arg->data,&result,1,MPI_INT, MPI_MIN,0, OP_MPI_WORLD);
	    memcpy(arg->data, &result, sizeof(int));
	}
    }
}



/**-------------------Performance measurement and reporting------------------**/

void op_mpi_timing_output(int my_rank)
{
    printf("\n\nPerformance information on rank %d\n", my_rank);
    
    if (OP_kern_max>0) {
    printf("\n  count     time     kernel name ");
    printf("\n -------------------------------- \n");
    for (int n=0; n<OP_kern_max; n++) {
      if (OP_kernels[n].count>0) {
          printf(" %6d  %8.4f      %s \n",
	       OP_kernels[n].count,
               OP_kernels[n].time,
               OP_kernels[n].name);
      }
    }
  }
}



/**-------------------------------Debug functions----------------------------**/

//print all lists for a given set 
/*void print_set_list(op_set set, int rank, int exec, int import)
{
    //exec import list
    
    //exec export list
    if(exec==1 & import ==0)
    {
    	set_halo_list l= OP_export_sets_list[set->index];
    	printf("On rank %d Exec export list for set %10s\n",rank, set->name);
    	for(int j=0; j<l->ranks_size;j++)
    	{
    	    for(int k=0;k<l->sizes[j];k++)
    	    printf("export to %d on rank %d set %6s contains %d\n",
    	    	l->ranks[j],rank,l->set->name,l->list[l->disps[j]+k]);
    	}
    }
    
    //nonexec import list
    
    
    //nonexec export list
    
}*/
