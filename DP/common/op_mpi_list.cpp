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



#include "op_mpi_datatypes.h"
#include "op_mpi_util.cpp"
#include "timers.h"

MPI_Comm OP_MPI_WORLD;

/*-----------------------MPI related global variables -------------------------*/
int OP_export_maps_index=0; 
int OP_export_sets_index=0;

int OP_import_maps_index=0; 
int OP_import_sets_index=0;

int OP_import_nonexec_set_index=0;
int OP_export_nonexec_set_index=0;


/*------------------------------MPI related global lists-----------------------------*/

map_export_list  * OP_export_maps_list[10]; 
set_export_list  * OP_export_sets_list[10]; 

map_import_list  * OP_import_maps_list[10]; 
set_import_list  * OP_import_sets_list[10];

nonexec_set_import_list  * OP_import_nonexec_sets_list[10];
nonexec_set_export_list  * OP_export_nonexec_sets_list[10];  


//declare lists  
map_export_list mxlist[10]; //map entries for export exec elements 
set_export_list sxlist[10]; //set export exec elements
map_import_list milist[10]; //map entries for import exec elements
set_import_list silist[10]; //set import exec elements
nonexec_set_export_list nesxlist[10]; //set export non-exec elements
nonexec_set_import_list neilist[10]; //set import non-exec elements
  

//global allray to hold dirty_bits for op_dats
int dirtybit[10];


//
//debug functions
//
#if DEBUG
#include "op_mpi_debug.cpp"
#endif

/*----------------------MPI list creation related functions-------------------*/

//declare export list for a given mapping table 
void decl_map_export_list(op_map map, int size, int* ranks, int ranks_size, 
	int* disps, int* sizes,	int* exp_list,  map_export_list &list){
  list.map = map;
  list.size = size;
  list.ranks = ranks;
  list.ranks_size = ranks_size;
  list.disps = disps;
  list.sizes = sizes;
  list.exp_list = exp_list;
  OP_export_maps_list[map.index] = &list;
  OP_export_maps_index++;
}

//declare export list for a set 
void decl_set_export_list(op_set set, int size, int* ranks, int ranks_size, 
	int* disps, int* sizes,	int* exp_list,  set_export_list &list){
  list.set = set;
  list.size = size;
  list.ranks = ranks;
  list.ranks_size = ranks_size;
  list.disps = disps;
  list.sizes = sizes;
  list.exp_list = exp_list;
  OP_export_sets_list[set.index] = &list;
  OP_export_sets_index++;
}

//declare import list for a given mapping table 
void decl_map_import_list(op_map map, int size, int* ranks, int ranks_size, 
	int* disps, int* sizes,	int* imp_list,  map_import_list &list){
  list.map = map;
  list.size = size;
  list.ranks = ranks;
  list.ranks_size = ranks_size;
  list.disps = disps;
  list.sizes = sizes;
  list.imp_list = imp_list;
  OP_import_maps_list[map.index] = &list;
  OP_import_maps_index++;
}

//declare import list for a set 
void decl_set_import_list(op_set set, int size, int* ranks, int ranks_size, 
	int* disps, int* sizes,	int* imp_list,  set_import_list &list){
  list.set = set;
  list.size = size;
  list.ranks = ranks;
  list.ranks_size = ranks_size;
  list.disps = disps;
  list.sizes = sizes;
  list.imp_list = imp_list;
  OP_import_sets_list[set.index] = &list;
  OP_import_sets_index++;
}

//declare nonexec import list for a set 
void decl_nonexec_set_import_list(op_set set, int size, int* ranks, int ranks_size, 
	int* disps, int* sizes,	int* imp_list,  nonexec_set_import_list &list){
  list.set = set;
  list.size = size;
  list.ranks = ranks;
  list.ranks_size = ranks_size;
  list.disps = disps;
  list.sizes = sizes;
  list.imp_list = imp_list;
  OP_import_nonexec_sets_list[set.index] = &list;
  OP_import_nonexec_set_index++;
  
}

//declare nonexec export list for a set 
void decl_nonexec_set_export_list(op_set set, int size, int* ranks, int ranks_size, 
	int* disps, int* sizes,	int* exp_list,  nonexec_set_export_list &list){
  list.set = set;
  list.size = size;
  list.ranks = ranks;
  list.ranks_size = ranks_size;
  list.disps = disps;
  list.sizes = sizes;
  list.exp_list = exp_list;
  OP_export_nonexec_sets_list[set.index] = &list;
  OP_export_nonexec_set_index++;
}

//sort list with [rank,element] into list format
void create_map_export_list(op_map map, int* temp_list, map_export_list &list, 
	int size, int comm_size, int my_rank)
{
    int* ranks = (int *) malloc(comm_size*sizeof(int));
    int* exp_list = (int *) malloc((size/2)*sizeof(int));
    int* disps = (int *) malloc(comm_size*sizeof(int));
    int* sizes = (int *) malloc(comm_size*sizeof(int));
    
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
	    //sort temp
	    quickSort(temp,0,sizes[index]-1);
	    
	    //eliminate duplicates in temp
	    sizes[index] = removeDups(temp, sizes[index]);
	    total_size = total_size + sizes[index];
	    
	    if(index > 0)
	    	disps[index] = disps[index-1] +  sizes[index-1];
	    //add to end of exp_list
	    for(int e = 0;e<sizes[index];e++)
	    	exp_list[disps[index]+e] = temp[e];
	    
	    index++;
	}
	free(temp);
    }
    decl_map_export_list(map, total_size, ranks, index, disps, sizes, exp_list, list);
}


void create_set_export_list(op_set set, int* temp_list, set_export_list &list, 
	int size, int comm_size, int my_rank)
{
    int* ranks = (int *) malloc(comm_size*sizeof(int));
    int* exp_list = (int *) malloc((size/2)*sizeof(int));
    int* disps = (int *) malloc(comm_size*sizeof(int));
    int* sizes = (int *) malloc(comm_size*sizeof(int));
    
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
    	    	exp_list[disps[index]+e] = temp[e];
    	    
    	    index++;
    	}
	free(temp);
    }
    decl_set_export_list(set, total_size, ranks, index, disps, sizes, exp_list, list);
}

void create_map_import_list(op_map map, int* temp_list, map_import_list &list, 
	int size, int* ranks, int* sizes, int ranks_size, int comm_size, int my_rank)
{
    int* disps = (int *) malloc(comm_size*sizeof(int));
    disps[0] = 0;
    for(int i=0; i<ranks_size; i++)
    {
    	if(i>0)disps[i] = disps[i-1]+sizes[i-1]; 	
    }
    decl_map_import_list(map, size, ranks, ranks_size, disps, sizes, temp_list,list);
}

void create_set_import_list(op_set set, int* temp_list, set_import_list &list, 
	int size, int* ranks, int* sizes, int ranks_size, int comm_size, int my_rank)
{
    int* disps = (int *) malloc(comm_size*sizeof(int));
    disps[0] = 0;
    for(int i=0; i<ranks_size; i++)
    {
    	if(i>0)disps[i] = disps[i-1]+sizes[i-1];
    }
    decl_set_import_list(set, size, ranks, ranks_size, disps, sizes, temp_list,list);
}


void create_nonexec_set_import_list(op_set set, int* temp_list, nonexec_set_import_list &list,
    int size, int comm_size, int my_rank)
{
    int* ranks = (int *) malloc(comm_size*sizeof(int));
    int* imp_list = (int *) malloc((size/2)*sizeof(int));
    int* disps = (int *) malloc(comm_size*sizeof(int));
    int* sizes = (int *) malloc(comm_size*sizeof(int));      	
    
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
	sizes[index] = 0;
	disps[index] = 0;
	int* temp = (int *) malloc((size/2)*sizeof(int)); 
	for(int i = 0;i<size;i=i+2)
	{
	    if(temp_list[i]==r)
	    	temp[sizes[index]++] = temp_list[i+1];	    
	}
	
	if(sizes[index]>0) 
	{
	    ranks[index] = r;
	    //sort temp
	    quickSort(temp,0,sizes[index]-1);
	    
	    //eliminate duplicates in temp
	    sizes[index] = removeDups(temp, sizes[index]);
	    total_size = total_size + sizes[index];
	    
	    if(index > 0) 
	    	disps[index] = disps[index-1] +  sizes[index-1];
	    
	    //add to end of exp_list
	    for(int e = 0;e<sizes[index];e++)
	    	imp_list[disps[index]+e] = temp[e];

	    index++;
		
	}
	free(temp);
    }
    decl_nonexec_set_import_list(set, total_size, ranks, index, disps, sizes, imp_list,list);
}

void create_nonexec_set_export_list(op_set set, int* temp_list, nonexec_set_export_list &list, 
	int size, int* ranks, int* sizes, int ranks_size, int comm_size, int my_rank)
{
    int* disps = (int *) malloc(comm_size*sizeof(int));
    disps[0] = 0;
    for(int i=0; i<ranks_size; i++)
    {
    	if(i>0)disps[i] = disps[i-1]+sizes[i-1];
    }
    decl_nonexec_set_export_list(set, size, ranks, ranks_size, disps, sizes, temp_list,list);
}

template <class T>
void find_neighbors(T* List, int* neighbors, int* sizes, 
	int* ranks_size, int my_rank, int comm_size)
{
    //map_export_list list=*OP_export_maps_list[map.index];
    int *temp = allocate(temp,comm_size*sizeof(int));
    int *r_temp = allocate(r_temp,comm_size*comm_size*sizeof(int));
    
    for(int r = 0;r<comm_size*comm_size;r++)r_temp[r] = -99;
    for(int r = 0;r<comm_size;r++)temp[r] = -99;    
    
    int n = 0;
    
    for(int r =0; r<comm_size; r++)
    {
    	if(List->ranks[r]>=0) temp[List->ranks[r]] = List->sizes[r];
    }
    
    MPI_Allgather( temp, comm_size, MPI_INT, r_temp,
    	comm_size,MPI_INT,OP_MPI_WORLD);

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

template <class T>
void get_neighbors(T* List, int* neighbors, int* sizes, 
	int* ranks_size, int my_rank, int comm_size)
{
    int *temp = allocate(temp,2*comm_size*sizeof(int));
    int n = 0;
    
    for(int r=0; r<comm_size; r++) {
    	if(my_rank==r)
    	{
    	    memcpy(&temp[0], (void *)&List->ranks[0],comm_size*sizeof(int));
    	    memcpy(&temp[comm_size], (void *)&List->sizes[0],comm_size*sizeof(int));
    	}
    	MPI_Bcast(temp, 2*comm_size, MPI_INT, r, OP_MPI_WORLD );
    	if(my_rank!=r)
    	{
    	    for(int i=0; i<comm_size; i++)
    	    {
    	    	if(temp[i] == my_rank)
    	    	{
    	    	    neighbors[n] = r;
    	    	    sizes[n] = temp[i+comm_size];
    	    	    n++;
    	    	}
    	    }
    	}
    }
    *ranks_size = n;
    free(temp);
}


 
/*--------------------Create Import Export Lists for MPI----------------------*/

void op_list_create()
{
  //declare timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;                                        
  double time;
  double max_time;
    
  timers(&cpu_t1, &wall_t1); //timer start for list create
  
  
  //create new communicator for OP mpi operation
  int my_rank, comm_size;
  MPI_Comm_dup(MPI_COMM_WORLD, &OP_MPI_WORLD);
  MPI_Comm_rank(OP_MPI_WORLD, &my_rank);
  MPI_Comm_size(OP_MPI_WORLD, &comm_size);
    
  
/**-------------------BEGIN Import/Export List Creation ---------------------**/

    /* Compute global partition range information for each set*/
    int** part_range = (int **) allocate(part_range,OP_set_index*sizeof(int*));
    get_part_range(part_range,comm_size);

  
    /*----- STEP 1 - construct export lists for mappings and execute sets-------*/
    
    //declare temporaty scratch variables to hold set export lists and mapping 
    //table export lists
    int s_i;
    int* set_list;
  
    int* map_list[OP_map_index];
    int cap_s = 1000; //keep track of the temp array capacities
  
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=*OP_set_list[s];
    	
    	//create a temporaty scratch space to hold export list for this set
    	s_i = 0;
    	set_list = allocate(set_list, cap_s*sizeof(int));
    	for(int e=0; e<set.size;e++){//for each elment of this set
    	    for(int m=0; m<OP_map_index; m++) { //for each maping table
    	    	op_map map=*OP_map_list[m];
    	    	
    	    	if(compare_sets(map.from,set)==1) //need to select mappings FROM this set
    	    	{
    	    	    int part, local_index;
    	    	    for(int j=0; j<map.dim; j++) { //for each element pointed at by this entry
    	    	    	part = get_partition(map.map[e*map.dim+j],
    	    	    	    part_range[map.to.index],&local_index,comm_size);
    	    	    	if(s_i>=cap_s)
    	    	    	{
    	    	    	    cap_s = cap_s*2;
    	    	    	    set_list = reallocate(set_list, cap_s*sizeof(int));
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
    	//printf("creating set export list for set %10s of size %d\n",set.name,s_i);
    	create_set_export_list(set,set_list, sxlist[set.index], s_i, comm_size, my_rank);
    	
    	for(int m=0; m<OP_map_index; m++) { //for each maping table create export list
    	    op_map map=*OP_map_list[m];
    	    if(compare_sets(map.from,set)==1) //need to select mappings FROM this set
    	    {
    	    	//create map export list: the union of all entries in mapping tables that is FROM this set
    	    	//printf("creating map export list for map %10s of size %d\n",map.name,s_i);
    	    	create_map_export_list(map,set_list, mxlist[map.index], s_i, comm_size, my_rank);
    	    	
    	    }
    	}
    	free(set_list);//free temp list
    }
    
       
    
    /*----- STEP 2 - construct import lists for mappings and execute sets-------*/
    
    /* Mappings first................*/
    
    int *neighbors, *sizes;
    int ranks_size;
    
    for(int m=0; m<OP_map_index; m++) { //for each maping table
    	op_map map=*OP_map_list[m];
    	
    	//-----Discover neighbors-----
    	ranks_size = 0;
    	neighbors = allocate(neighbors,comm_size*sizeof(int));
    	sizes = allocate(sizes,comm_size*sizeof(int));
    	
    	map_export_list list=*OP_export_maps_list[map.index];
    	
    	//find_neighbors(&list,neighbors,sizes,&ranks_size,my_rank,comm_size);
    	get_neighbors(&list,neighbors,sizes,&ranks_size,my_rank,comm_size);
    	MPI_Request request_send[list.ranks_size];
    	
    	int* rbuf, cap = 0, index = 0;
    	
    	for(int i=0; i< ranks_size; i++) cap = cap + sizes[i];    	     
    	int* temp = allocate(temp,cap*sizeof(int));
 
    	for(int i=0; i<list.ranks_size; i++) {
    	    //printf("export map %10s to %d from rank %d, list of size %d \n",
    	    //    map.name,list.ranks[i],my_rank,list.sizes[i]);
    	    int* sbuf = &list.exp_list[list.disps[i]];
    	    MPI_Isend( sbuf,  list.sizes[i],  MPI_INT, list.ranks[i], m,
    	    	OP_MPI_WORLD, &request_send[i]);
    	}
     
    	//import this list from those neighbors
    	for(int i=0; i<ranks_size; i++) {
    	    //printf("on rank %d, import map %10s from %d list of size %d \n",
    	    //    my_rank,map.name,neighbors[i],sizes[i]);
    	    rbuf = allocate(rbuf,sizes[i]*sizeof(int));
    	    MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], m, OP_MPI_WORLD,
    	    	MPI_STATUSES_IGNORE );
    	    memcpy(&temp[index],(void *)&rbuf[0],sizes[i]*sizeof(int));
    	    index = index + sizes[i];
    	    free(rbuf);
    	}
 
    	MPI_Waitall(list.ranks_size,request_send, MPI_STATUSES_IGNORE );
    	
    	//create import lists
    	//printf("creating importlist of with number of neighbors %d\n",ranks_size);
    	create_map_import_list(map, temp, milist[map.index],index,neighbors, 
    	    sizes, ranks_size, comm_size, my_rank);
    }
    
        
    /* Sets next..............*/
    
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=*OP_set_list[s];
    	
    	//-----Discover neighbors-----
    	ranks_size = 0;
    	neighbors = allocate(neighbors,comm_size*sizeof(int));
    	sizes = allocate(sizes,comm_size*sizeof(int));
    	
    	set_export_list list=*OP_export_sets_list[set.index];
    	
    	//find_neighbors(&list,neighbors,sizes,&ranks_size,my_rank,comm_size);
    	get_neighbors(&list,neighbors,sizes,&ranks_size,my_rank,comm_size);
    	
    	MPI_Request request_send[list.ranks_size];
    	
    	int* rbuf, cap = 0, index = 0;
    	
    	for(int i=0; i<list.ranks_size; i++) {
    	    //printf("export from %d to %d set %10s, list of size %d \n",
    	    //my_rank,list.ranks[i],set.name,list.sizes[i]);
    	    int* sbuf = &list.exp_list[list.disps[i]];
    	    MPI_Isend( sbuf,  list.sizes[i],  MPI_INT, list.ranks[i], s,
    	    	OP_MPI_WORLD, &request_send[i] );
    	}
    	
    	for(int i=0; i< ranks_size; i++) cap = cap + sizes[i];
    	int* temp = allocate(temp,cap*sizeof(int));
    	
    	//import this list from those neighbors
    	for(int i=0; i<ranks_size; i++) {
    	    //printf("import from %d to %d set %10s, list of size %d\n",
    	    //neighbors[i], my_rank, set.name, sizes[i]);
    	    rbuf = allocate(rbuf,sizes[i]*sizeof(int));
    	    MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i],s, OP_MPI_WORLD,
    	    	MPI_STATUSES_IGNORE );
    	    memcpy(&temp[index],(void *)&rbuf[0],sizes[i]*sizeof(int));
    	    index = index + sizes[i];
    	    free(rbuf);
    	}
    	
    	MPI_Waitall(list.ranks_size,request_send, MPI_STATUSES_IGNORE );
    	
    	//create import lists
    	//printf("creating importlist with number of neighbors %d\n",ranks_size);
    	create_set_import_list(set, temp, silist[set.index],index,neighbors,
    	    sizes, ranks_size, comm_size, my_rank);
    }
  
    
    
 
    /*-- STEP 3 - Exchange mapping table entries using the import/export lists--*/
  
    for(int m=0; m<OP_map_index; m++) { //for each maping table
    	op_map map=*OP_map_list[m];
    	map_import_list i_list = *OP_import_maps_list[map.index];
    	map_export_list e_list = *OP_export_maps_list[map.index];
  	  
    	MPI_Request request_send[e_list.ranks_size];

    	//prepare bits of the mapping tables to be exported
    	for(int i=0; i<e_list.ranks_size; i++) {
    	    int* sbuf = allocate(sbuf,e_list.sizes[i]*map.dim*sizeof(int));
    	    for(int j = 0; j<e_list.sizes[i]; j++)
    	    {
    	    	for(int p = 0; p< map.dim; p++)
    	    	{
    	    	    sbuf[j*map.dim+p] =
    	    	    map.map[map.dim*(e_list.exp_list[e_list.disps[i]+j])+p];
    	    	}
    	    }
    	    //printf("\n export from %d to %d map %10s, number of elements of size %d | sending:\n ",
    	    //    my_rank,e_list.ranks[i],map.name,e_list.sizes[i]);
    	    MPI_Isend(sbuf,  map.dim*e_list.sizes[i],  MPI_INT, e_list.ranks[i],
    	    	m, OP_MPI_WORLD, &request_send[i]);
    	}
      
    	//prepare space for the incomming mapping tables - realloc each
    	//mapping tables in each mpi process
    	OP_map_list[map.index]->map = reallocate(OP_map_list[map.index]->map,
    	    (map.dim*(map.from.size+i_list.size))*sizeof(int));
      
    	int init = map.dim*(map.from.size);
    	for(int i=0; i<i_list.ranks_size; i++) {
    	    MPI_Recv(&(OP_map_list[map.index]->map[init+i_list.disps[i]*map.dim]),
    	    	map.dim*i_list.sizes[i], MPI_INT, i_list.ranks[i], m,
    	    	OP_MPI_WORLD, MPI_STATUSES_IGNORE);
    	}

    	MPI_Waitall(e_list.ranks_size,request_send, MPI_STATUSES_IGNORE );
    	//printf("\n imported on to %d map %10s, number of elements of size %d | recieving: ",
    	//	  my_rank, map.name, i_list.size);
    }
    
    
    /*-- STEP 4 - Create import lists for non-execute set elements using mapping
    table entries including the additional mapping table entries --*/
    
    //declare temporaty scratch variables to hold non-exec set export lists
    s_i = 0;
    set_list = NULL;  
    cap_s = 1000; //keep track of the temp array capacity
  
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=*OP_set_list[s];
    	set_import_list exec_set_list=*OP_import_sets_list[set.index];   
      
    	//create a temporaty scratch space to hold nonexec export list for this set
    	s_i = 0;
    	set_list = allocate(set_list, cap_s*sizeof(int));     
      
    	for(int m=0; m<OP_map_index; m++) { //for each maping table
    	    op_map map=*OP_map_list[m];
    	    map_import_list exec_map_list=*OP_import_maps_list[map.index];
 
    	    if(compare_sets(map.to,set)==1) //need to select mappings TO this set
    	    {
    	    	//for each entry in this mapping table: original+execlist
    	    	int len = map.from.size+exec_map_list.size;
    	    	for(int e = 0; e<len; e++)
    	    	{
    	    	    int part;
    	    	    int local_index;
    	    	    for(int j=0; j<map.dim; j++) { //for each element pointed at by this entry
    	    	    	part = get_partition(map.map[e*map.dim+j],
    	    	    	    part_range[map.to.index],&local_index,comm_size);
              	      
    	    	    	if(s_i>=cap_s)
    	    	    	{
    	    	    	    cap_s = cap_s*2;
    	    	    	    set_list = reallocate(set_list,cap_s*sizeof(int));
    	    	    	}
              	      
    	    	    	if(part != my_rank)
    	    	    	{
    	    	    	    int found = -1;
    	    	    	    //check in exec list
    	    	    	    int rank = binary_search(exec_set_list.ranks,
    	    	    	    	part, 0, exec_set_list.ranks_size-1);
              	      	  
    	    	    	    if(rank >= 0)
    	    	    	    {
    	    	    	    	found = binary_search(exec_set_list.imp_list,
    	    	    	    	    local_index, exec_set_list.disps[rank],
    	    	    	    	    exec_set_list.disps[rank]+exec_set_list.sizes[rank]-1);
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
    	create_nonexec_set_import_list(set,set_list, neilist[set.index], s_i, comm_size, my_rank);
    	free(set_list);//free temp list
      
    }
  
  
    /*----------- STEP 5 - construct non-execute set export lists -------------*/
  
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=*OP_set_list[s];
      
    	//-----Discover neighbors-----
    	ranks_size = 0;
    	neighbors = allocate(neighbors,comm_size*sizeof(int));
    	sizes = allocate(sizes,comm_size*sizeof(int));
      
    	nonexec_set_import_list list=*OP_import_nonexec_sets_list[set.index];
    	get_neighbors(&list,neighbors,sizes,&ranks_size,my_rank,comm_size);
    	//find_neighbors(&list,neighbors,sizes,&ranks_size,my_rank,comm_size);
            
    	MPI_Request request_send[list.ranks_size];
            
    	int* rbuf, cap = 0, index = 0;
      
    	for(int i=0; i<list.ranks_size; i++) {
    	    //printf("import to %d from %d set %10s, nonexec list of size %d | sending:\n",
    	    //    my_rank,list.ranks[i],set.name,list.sizes[i]);
    	    int* sbuf = &list.imp_list[list.disps[i]];
    	    MPI_Isend( sbuf,  list.sizes[i],  MPI_INT, list.ranks[i], s,
    	    	OP_MPI_WORLD, &request_send[i] );  
    	}
      
    	for(int i=0; i< ranks_size; i++) cap = cap + sizes[i];     	     
    	int* temp = allocate(temp,cap*sizeof(int));
     
    	//export this list to those neighbors
    	for(int i=0; i<ranks_size; i++) {
    	    //printf("export to %d from %d set %10s, list of size %d | recieving:\n",
    	    //    neighbors[i], my_rank, set.name, sizes[i]);
    	    rbuf = allocate(rbuf,sizes[i]*sizeof(int));
    	    MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i],s, OP_MPI_WORLD,
    	    	MPI_STATUSES_IGNORE );
    	    memcpy(&temp[index],(void *)&rbuf[0],sizes[i]*sizeof(int));
    	    index = index + sizes[i];
    	    free(rbuf);
    	}
      
    	MPI_Waitall(list.ranks_size,request_send, MPI_STATUSES_IGNORE );
     
    	//create import lists
    	//printf("creating nonexec set export list with number of neighbors %d\n",ranks_size);
    	create_nonexec_set_export_list(set, temp, nesxlist[set.index],index,
    	    neighbors, sizes, ranks_size, comm_size, my_rank);
    }
    
/**-------------------- END Import/Export List Creation ---------------------**/
  


  
/**-BEGIN Initial Exchange of All Actual Halo Data Using Import/Export Lists-**/

    /*-STEP 1 - Exchange execute set elements/data using the import/export lists--*/

    for(int s=0; s<OP_set_index; s++){ //for each set
    	op_set set=*OP_set_list[s];
    	set_import_list i_list = *OP_import_sets_list[set.index];
    	set_export_list e_list = *OP_export_sets_list[set.index];
      
    	//for each data array
    	for(int d=0; d<OP_dat_index; d++){
    	    op_dat dat=*OP_dat_list[d];
      	  
    	    if(compare_sets(set,dat.set)==1)//if this data array is defined on this set
    	    {
    	    	//printf("on rank %d, The data array is %10s\n",my_rank,dat.name);
    	    	MPI_Request request_send[e_list.ranks_size];
      	      
    	    	//prepare execute set element data to be exported
    	    	for(int i=0; i<e_list.ranks_size; i++) {
    	    	    char* sbuf = allocate(sbuf,e_list.sizes[i]*dat.size);
    	    	    for(int j = 0; j<e_list.sizes[i]; j++)
    	    	    {
    	    	    	int set_elem_index = e_list.exp_list[e_list.disps[i]+j];
    	    	    	memcpy(&sbuf[j*dat.size],(void *)&dat.dat[dat.size*(set_elem_index)],dat.size);
    	    	    	/*for(int v = 0; v< dat.size; v++)
    	    	    	{
    	    	    	    sbuf[j*dat.size+v] = dat.dat[dat.size*(set_elem_index)+v];
    	    	    	}*/
    	    	    }
    	    	    //printf("export from %d to %d data %10s, number of elements of size %d | sending:\n ",
    	    	    //    my_rank,e_list.ranks[i],dat.name,e_list.sizes[i]);
    	    	    MPI_Isend(sbuf,  dat.size*e_list.sizes[i],  MPI_CHAR, e_list.ranks[i],
    	    	    	d, OP_MPI_WORLD, &request_send[i]);
    	    	}
      	      
    	    	//prepare space for the incomming data - realloc each
    	    	//data array in each mpi process
    	    	OP_dat_list[dat.index]->dat = reallocate(OP_dat_list[dat.index]->dat,
    	    	    (set.size+i_list.size)*dat.size);
      	      
    	    	int init = set.size*dat.size;
      	      
    	    	for(int i=0; i<i_list.ranks_size; i++) {
    	    	    MPI_Recv(&(OP_dat_list[dat.index]->dat[init+i_list.disps[i]*dat.size]),
    	    	    	dat.size*i_list.sizes[i], MPI_CHAR, i_list.ranks[i], d,
    	    	    	OP_MPI_WORLD, MPI_STATUSES_IGNORE);
    	    	}
      	      
    	    	MPI_Waitall(e_list.ranks_size,request_send, MPI_STATUSES_IGNORE );
    	    	//printf("imported on to %d data %10s, number of elements of size %d | recieving:\n ",
    	    	//	  my_rank, dat.name, i_list.size);
    	    }
    	}
    }


    /*-STEP 2 - Exchange non-execute set elements/data using the import/export lists--*/

    for(int s=0; s<OP_set_index; s++){ //for each set
    	op_set set=*OP_set_list[s];
    	nonexec_set_import_list i_list = *OP_import_nonexec_sets_list[set.index];
    	nonexec_set_export_list e_list = *OP_export_nonexec_sets_list[set.index];
        
    	//for each data array
    	for(int d=0; d<OP_dat_index; d++){
    	    op_dat dat=*OP_dat_list[d];
      	  
    	    if(compare_sets(set,dat.set)==1)//if this data array is defined on this set
    	    {
    	    	//printf("on rank %d, The data array is %10s\n",my_rank,dat.name);
    	    	MPI_Request request_send[e_list.ranks_size];
      	      
    	    	//prepare non-execute set element data to be exported
    	    	for(int i=0; i<e_list.ranks_size; i++) {
    	    	    char* sbuf = allocate(sbuf,e_list.sizes[i]*dat.size);
    	    	    for(int j = 0; j<e_list.sizes[i]; j++)
    	    	    {
    	    	    	int set_elem_index = e_list.exp_list[e_list.disps[i]+j];
    	    	    	memcpy(&sbuf[j*dat.size],(void *)&dat.dat[dat.size*(set_elem_index)],dat.size);
    	    	    	/*for(int v = 0; v< dat.size; v++)
    	    	    	{
    	    	    	    sbuf[j*dat.size+v] = dat.dat[dat.size*(set_elem_index)+v];
    	    	    	}*/
    	    	    }
    	    	    MPI_Isend(sbuf,  dat.size*e_list.sizes[i],  MPI_CHAR, e_list.ranks[i],
    	    	    	d, OP_MPI_WORLD, &request_send[i]);
    	    	}
      	      
    	    	//prepare space for the incomming nonexec-data - realloc each
    	    	//data array in each mpi process
      	      
    	    	set_import_list exec_i_list = *OP_import_sets_list[set.index];      
      
    	    	OP_dat_list[dat.index]->dat = reallocate(OP_dat_list[dat.index]->dat,
    	    	    (set.size+exec_i_list.size+i_list.size)*dat.size);
    	    	int init = (set.size+exec_i_list.size)*dat.size;
      	      
    	    	for(int i=0; i<i_list.ranks_size; i++) {
    	    	    MPI_Recv(&(OP_dat_list[dat.index]->dat[init+i_list.disps[i]*dat.size]),
    	    	    	dat.size*i_list.sizes[i], MPI_CHAR, i_list.ranks[i], d,
    	    	    	OP_MPI_WORLD, MPI_STATUSES_IGNORE);
    	    	}
      	      
    	    	MPI_Waitall(e_list.ranks_size,request_send, MPI_STATUSES_IGNORE );
    	    }
    	}
    }

  
  
    /*-STEP 3 ------------------- Renumber Mapping tables-------------------------*/
	
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=*OP_set_list[s];
      
    	for(int m=0; m<OP_map_index; m++) { //for each maping table
    	    op_map map=*OP_map_list[m];
      	        	 
    	    if(compare_sets(map.to,set)==1) //need to select mappings TO this set
    	    {	
    	    	set_import_list exec_set_list=*OP_import_sets_list[set.index];
    	    	nonexec_set_import_list nonexec_set_list=*OP_import_nonexec_sets_list[set.index];
      	  
    	    	map_import_list exec_map_list=*OP_import_maps_list[map.index];
              
    	    	//for each entry in this mapping table: original+execlist
    	    	int len = map.from.size+exec_map_list.size;
    	    	for(int e = 0; e<len; e++)
    	    	{
    	    	    for(int j=0; j<map.dim; j++) { //for each element pointed at by this entry
    	    	    	
    	    	    	int part;
    	    	    	int local_index = 0;
    	    	    	part = get_partition(map.map[e*map.dim+j],
    	    	    	    part_range[map.to.index],&local_index,comm_size);
              	                    	      
    	    	    	if(part == my_rank)
    	    	    	{
    	    	    	    OP_map_list[map.index]->map[e*map.dim+j] = local_index;
    	    	    	}
    	    	    	else
    	    	    	{
    	    	    	    int found = -1;
    	    	    	    //check in exec list
    	    	    	    int rank1 = binary_search(exec_set_list.ranks,
    	    	    	    	part, 0, exec_set_list.ranks_size-1);
    	    	    	    //check in nonexec list
    	    	    	    int rank2 = binary_search(nonexec_set_list.ranks,
    	    	    	    	part, 0, nonexec_set_list.ranks_size-1);
              	      	  
    	    	    	    if(rank1 >=0)
    	    	    	    {
    	    	    	    	found = binary_search(exec_set_list.imp_list,
    	    	    	    	    local_index, exec_set_list.disps[rank1],
    	    	    	    	    exec_set_list.disps[rank1]+exec_set_list.sizes[rank1]-1);

              	      	      	if(found>=0) 
              	      	      	{
              	      	      	    OP_map_list[map.index]->map[e*map.dim+j] = 
              	      	      	    found + map.to.size ;
              	      	      	}
              	      	    }
              	      	  
              	      	    if(rank2 >=0 & found <0)
              	      	    {
              	      	    	found = binary_search(nonexec_set_list.imp_list,
              	      	    	    local_index, nonexec_set_list.disps[rank2],
              	      	    	    nonexec_set_list.disps[rank2]+
              	      	    	    nonexec_set_list.sizes[rank2]-1);
              	      	    	
              	      	    	if(found>=0)
              	      	      	{
              	      	      	    OP_map_list[map.index]->map[e*map.dim+j] = 
              	      	      	    found + set.size + exec_set_list.size;
              	      	      	}
              	      	    }
              	      	  
              	      	    if(found < 0) 
              	      	    printf("ERROR: Set %10s Element %d needed on rank %d \
              	      	    	from partition %d\n",set.name, local_index, my_rank,part );
              	      	}
              	    }
              	}
            }
        }
    }
    
    //set dirty bits of all data arrays to 0
    //for each data array
    for(int d=0; d<OP_dat_index; d++){
    	op_dat dat=*OP_dat_list[d];
    	
    	dirtybit[dat.index] = 0;
    }
    
    timers(&cpu_t2, &wall_t2);  //timer stop for list create
    
    //printf time to create import/export lists
    time = wall_t2-wall_t1;
    MPI_Reduce(&time,&max_time,1,MPI_DOUBLE, MPI_MAX,0, OP_MPI_WORLD);
    if(my_rank==0)printf("Max total halo creation time = %lf\n",max_time);
}






