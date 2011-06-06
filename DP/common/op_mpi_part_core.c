//parmetis header
#include <parmetis.h>

#include "op_mpi_part_core.h"

//
//MPI Communicator for partitioning
//
MPI_Comm OP_PART_WORLD;


/**----------------MPI partitioning related global variables ----------------**/
part* OP_part_list;
part* OP_part_list_original; 


/**----------------MPI partitioning related utility functions ---------------**/

//declare partition information for a given set 
void decl_partition(op_set set, int* g_index, int* partition)
{
    
  part p = (part) malloc(sizeof(part_core));
  p->set = set;
  p->g_index = g_index;
  p->elem_part = partition;
  p->is_partitioned = 0;
  OP_part_list[set->index] = p;
  
}

int frequencyof(int value, int* array, int size)
{
   int frequency = 0; 
   for(int i = 0; i<size; i++)
   {
    	if(array[i] == value) frequency++;   
   }
   return frequency;
}

int find_mode(int* array, int size)
{
   int count = 0, mode = array[0], current;
   for(int i=0; i<size; i++)
   {
        current = frequencyof(array[i], array, size);
    	if(count< current)
    	{   	    	
    	    count = current;
    	    mode = array[i];
    	}
   }
   return mode;
}

int compare_all_sets(op_set target_set, op_set other_sets[], int size)
{
    for(int i = 0; i < size; i++)
    {
    	//printf("comparing set %s with set %s\n",target_set.name, other_sets[i].name);
 	if(compare_sets(target_set, other_sets[i])==1)return i;
    }    
    return -1;
}


void create_exp_list(op_set set, int* temp_list, set_halo_list halo_list, 
    int size, int comm_size, int my_rank)
{
    int* ranks = (int *) malloc(comm_size*sizeof(int));
    int* list = (int *) malloc((size/2)*sizeof(int));
    int* disps = (int *) malloc(comm_size*sizeof(int));
    int* sizes = (int *) malloc(comm_size*sizeof(int));
    
    int ranks_size = 0;
    int total_size = 0;  
    
    create_list(list, ranks, disps, sizes, &ranks_size, &total_size,
    temp_list, size, comm_size, my_rank);
    
    halo_list->set = set;
    halo_list->size = total_size;
    halo_list->ranks = ranks;
    halo_list->ranks_size = ranks_size;
    halo_list->disps = disps;
    halo_list->sizes = sizes;
    halo_list->list = list;
}

void create_imp_list(op_set set, int* temp_list, set_halo_list halo_list, int total_size, 
    	int* ranks, int* sizes, int ranks_size, int comm_size, int my_rank)
{
    int* disps = (int *) malloc(comm_size*sizeof(int));
    disps[0] = 0;
    for(int i=0; i<ranks_size; i++)
    {
    	if(i>0)disps[i] = disps[i-1]+sizes[i-1]; 	
    }

    halo_list->set = set;
    halo_list->size = total_size;
    halo_list->ranks = ranks;
    halo_list->ranks_size = ranks_size;
    halo_list->disps = disps;
    halo_list->sizes = sizes;
    halo_list->list = temp_list;
    
}


/** ----use the partitioned map->to set to partition the map->from set ------**/
int partition_from_set(op_map map, int my_rank, int comm_size, int** part_range)
{
    
    part p_set = OP_part_list[map->to->index];
    
    int cap = 100; int count = 0;
    int* temp_list = malloc(cap*sizeof(int));
    if(temp_list == NULL) {
    	printf(" partition_from_set -- error allocating memory: int* temp_list\n");
    	exit(-1);
    }
    
    set_halo_list pi_list = (set_halo_list) malloc(sizeof(set_halo_list_core));
    
    //go through the map and build an import list of the non-local "to" elements 
    for(int i = 0; i < map->from->size; i++)
    {
    	int part, local_index;
    	for(int j = 0; j<map->dim; j++)
    	{
    	    part = get_partition(map->map[i*map->dim+j],
    	    	    	    part_range[map->to->index],&local_index,comm_size);    	    
    	    if(count>=cap)
    	    {
    	    	cap = cap*2;
    	    	temp_list = realloc(temp_list,cap*sizeof(int));
    	    	if(temp_list == NULL) {
    	    	    printf(" partition_from_set -- error reallocating memory: int* temp_list\n");
    	    	    exit(-1);
    	    	}
    	    }
    	    
    	    if(part != my_rank)
    	    {
    	    	temp_list[count++] = part;
    	    	temp_list[count++] = local_index;
    	    }    	    
    	}    	
    }    
    create_exp_list(map->to,temp_list, pi_list, count, comm_size, my_rank);
    free(temp_list);
    
    
    //now, discover neighbors and create export list of "to" elements
    int ranks_size = 0;

    int* neighbors = malloc(comm_size*sizeof(int));
    if(neighbors == NULL) {
    	printf(" partition_from_set -- error allocating memory: int *neighbors\n");
    	exit(-1);
    }
    int* sizes = malloc(comm_size*sizeof(int));
    if(sizes == NULL) {
    	printf(" partition_from_set -- error allocating memory: int *sizes\n");
    	exit(-1);
    }
    
    set_halo_list pe_list = (set_halo_list) malloc(sizeof(set_halo_list_core));
    find_neighbors_set(pi_list,neighbors,sizes,&ranks_size,my_rank,comm_size, OP_PART_WORLD);
    
    MPI_Request request_send[pi_list->ranks_size];
    int* rbuf; 
    cap = 0; count = 0;
    
    for(int i=0; i < pi_list->ranks_size; i++) {
    	int* sbuf = &pi_list->list[pi_list->disps[i]];
    	MPI_Isend( sbuf,  pi_list->sizes[i],  MPI_INT, pi_list->ranks[i], 1,
    	    OP_PART_WORLD, &request_send[i] );
    }
    
    for(int i=0; i< ranks_size; i++) cap = cap + sizes[i];
   
    temp_list = malloc(cap*sizeof(int));
    if(temp_list == NULL) {
    	printf(" partition_from_set -- error allocating memory: int* temp_list\n");
    	exit(-1);
    }
    
    for(int i=0; i<ranks_size; i++) {
    	rbuf = malloc(sizes[i]*sizeof(int));
    	if(rbuf == NULL) {
    	    printf(" partition_from_set -- error allocating memory: int* rbuf\n");
    	    exit(-1);
    	}
    	    
    	MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], 1, OP_PART_WORLD,
    	    MPI_STATUSES_IGNORE );
    	memcpy(&temp_list[count],(void *)&rbuf[0],sizes[i]*sizeof(int));
    	count = count + sizes[i];
    	free(rbuf);
    }
    MPI_Waitall(pi_list->ranks_size,request_send, MPI_STATUSES_IGNORE );
    create_imp_list(map->to,temp_list, pe_list, count, 
    	neighbors, sizes, ranks_size, comm_size, my_rank);    
    
    
    //use the import and export lists to exchange partition information of 
    //this "to" set
    MPI_Request request_send_p[pe_list->ranks_size];
    
    //first - prepare partition information of the "to" set element to be exported
    int** sbuf = (int **)malloc(pe_list->ranks_size*sizeof(int *));
    if(sbuf == NULL) {
    	printf(" partition_from_set -- error allocating memory: int** sbuf\n");
    	exit(-1);
    }
    	
    for(int i = 0; i < pe_list->ranks_size; i++)
    {
    	//printf("export to %d from rank %d set %s of size %d\n",
    	//   pe_list->ranks[i], my_rank, map->to->name, pe_list->sizes[i] );
    	sbuf[i] = (int *)malloc(pe_list->sizes[i]*sizeof(int));
    	if(sbuf[i] == NULL) {
    	    printf(" partition_from_set -- error allocating memory: sbuf[i]\n");
    	    exit(-1);
    	}
    
    	for(int j = 0; j<pe_list->sizes[i]; j++)
    	{
    	    int elem = pe_list->list[pe_list->disps[i]+j];
    	    sbuf[i][j] = p_set->elem_part[elem];    	 	   
    	}
    	MPI_Isend(sbuf[i],  pe_list->sizes[i],  MPI_INT, pe_list->ranks[i],
    	    	    	2, OP_PART_WORLD, &request_send_p[i]);    	
    }
    
    //second - prepare space for the incomming partition information of the "to" set
    int* imp_part = (int *)malloc(sizeof(int)*pi_list->size); 
    if(imp_part == NULL) {
    	printf(" partition_from_set -- error allocating memory: int* imp_part\n");
    	exit(-1);
    }
    	
    //third - receive 
    for(int i=0; i<pi_list->ranks_size; i++) {
    	//printf("import from %d to rank %d set %s of size %d\n",
    	//    pi_list->ranks[i], my_rank, map->to->name, pi_list->sizes[i] );
    	MPI_Recv(&imp_part[pi_list->disps[i]],
    	    pi_list->sizes[i], MPI_INT, pi_list->ranks[i], 2, 
    	    OP_PART_WORLD, MPI_STATUSES_IGNORE);
    	
    }
    MPI_Waitall(pe_list->ranks_size,request_send_p, MPI_STATUSES_IGNORE );
    for(int i = 0; i<pe_list->ranks_size; i++) free(sbuf[i]);free(sbuf);
    
    //allocate memory to hold the partition details for the set thats going to be 
    //partitioned 
    int* partition = (int *)malloc(sizeof(int)*map->from->size);
    if(partition == NULL) {
    	printf(" partition_from_set -- error allocating memory: int* partition\n");
    	exit(-1);
    }
   
    //go through the mapping table and the imported partition information and
    //partition the "from" set
    for(int i = 0; i<map->from->size; i++)
    {
    	int part, local_index;
    	int found_parts[map->dim];
    	for(int j = 0; j < map->dim; j++)
    	{
    	    part = get_partition(map->map[i*map->dim+j],
    	    	    	    part_range[map->to->index],&local_index,comm_size);
    	    
    	    if(part == my_rank)
    	    	found_parts[j] = p_set->elem_part[local_index];
    	    else //get partition information from imported data  
    	    {
    	    	 int r = binary_search(pi_list->ranks,part,0,pi_list->ranks_size-1);
    	    	 if(r >= 0)
    	    	 {
    	    	     int elem = binary_search(&pi_list->list[pi_list->disps[r]],
    	    	     	 local_index,0,pi_list->sizes[r]-1);
    	    	     if(elem >= 0)
    	    	     	 found_parts[j] = imp_part[elem];
    	    	     else
    	    	     {
    	    	     	 printf("Element %d not found in partition import list\n",local_index);
    	    	     	 exit(1);
    	    	     }
    	    	 }
    	    	 else 
    	    	 {
    	    	     printf("Rank %d not found in partition import list\n", part);
    	    	     exit(1);
    	    	 }
    	    }    	 	   
    	}
	partition[i] = find_mode(found_parts, map->dim);     	
    }
            
    
    
    
    OP_part_list[map->from->index]->elem_part = partition;
    OP_part_list[map->from->index]->is_partitioned = 1;
    
    //cleanup
    free(imp_part);
    free(pi_list->list);free(pi_list->ranks);free(pi_list->sizes);
    free(pi_list->disps);free(pi_list);
    free(pe_list->list);free(pe_list->ranks);free(pe_list->sizes);
    free(pe_list->disps);free(pe_list);
    
    return 1;    
}


/** ---- use the partitioned map->from set to partition the map->to set -------**/
int partition_to_set(op_map map, int my_rank, int comm_size, int** part_range)
{
    //OP_part_info_list[map->to->index]->part = partition;
    OP_part_list[map->to->index]->is_partitioned = 1;
    return 1;
}



/**--------------Partition all secondary sets and data migration-------------**/   
void partition_all(op_set primary_set, int my_rank, int comm_size)
{
    
/*--STEP 1 - Partition Secondary sets using primary set partition ------------*/
    
    // Compute global partition range information for each set
    int** part_range = (int **)malloc(OP_set_index*sizeof(int*));
    if(part_range == NULL) {
    	printf(" partition_all -- error allocating memory: int** part_range\n");
    	exit(-1);
    }
    get_part_range(part_range,my_rank,comm_size, OP_PART_WORLD);

    int sets_partitioned = 1;
    int maps_used = 0;
    
    op_set all_partitioned_sets[OP_set_index]; 
    int all_used_maps[OP_map_index];
    for(int i = 0; i<OP_map_index; i++) { all_used_maps[i] = -1;}
    
    //begin with the partitioned primary set - e.g nodes
    all_partitioned_sets[0] = OP_set_list[primary_set->index];
    
    int error = 0;
    while(sets_partitioned < OP_set_index && error == 0)
    {
    	int cost[OP_map_index];
    	for(int i = 0; i<OP_map_index; i++) cost[i] = 99;
       
    	//compute a "cost" associated with using each mapping table 
    	for(int m=0; m<OP_map_index; m++)
    	{
    	    op_map map=OP_map_list[m];
      	   
    	    if(binary_search(all_used_maps,map->index,0,maps_used-1)<0)// if not used before 
    	    {
    	    	part to_set = OP_part_list[map->to->index];
    	    	part from_set = OP_part_list[map->from->index];
       	       
    	    	//partitioning a set using a mapping from a partitioned set costs
    	    	//more than partitioning a set using a mapping to a partitioned set
    	    	//i.e. preferance is given to the latter over the former       	       
    	    	if(from_set->is_partitioned == 1 && 
    	    	    compare_all_sets(map->from,all_partitioned_sets,sets_partitioned)>=0)
    	    	    cost[map->index] = 2;
    	    	else if(to_set->is_partitioned == 1 && 
    	    	    compare_all_sets(map->to,all_partitioned_sets,sets_partitioned)>=0)
    	    	    cost[map->index] = 0;
    	    }
       	}
       
       	while(1)
       	{
       	    int selected = min(cost, OP_map_index);
       	    if(selected >= 0)
       	    {
       	    	op_map map=OP_map_list[selected];
       	       
       	    	//partition using this map       	              	   
       	    	part to_set = OP_part_list[map->to->index];
       	    	part from_set = OP_part_list[map->from->index];
       	       
       	    	if(to_set->is_partitioned == 1) 
       	    	{
       	    	    if( partition_from_set(map, my_rank, comm_size, part_range) > 0)
       	    	    {
       	    	    	//if(my_rank==0)printf("Found map %s to partition from set %s using set %s\n",
       	    	    	//	   map->name,map->from->name,map->to->name);
       	    	    	all_partitioned_sets[sets_partitioned++] = map->from;
       	    	    	all_used_maps[maps_used++] = map->index;
       	    	    	break;       	       	       
       	    	    }
       	    	    else //partitioning unsuccessful with this map- find another map
       	    	    	cost[selected] = 99;
       	    	} 
       	    	else if(from_set->is_partitioned == 1) 
       	    	{
       	    	    if( partition_to_set(map, my_rank, comm_size, part_range) > 0)
       	    	    {
       	    	    	//if(my_rank==0)printf("Found map %s to partition to set %s using set %s\n",
       	    	    	//	   map->name,map->to->name,map->from->name);
       	    	    	all_partitioned_sets[sets_partitioned++] = map->to;
       	    	    	all_used_maps[maps_used++] = map->index;
       	    	    	break;
       	    	    }
       	    	    else //partitioning unsuccessful with this map - find another map
       	    	    	cost[selected] = 99;
       	    	}       	       
       	    }
       	    else //partitioning error;
       	    {
       	    	error = 1; break;
       	    }       	   
       	}
    }
       
    if(my_rank==0)
    {
    	printf("Sets partitioned = %d\n",sets_partitioned);   
    	if(sets_partitioned != OP_set_index)
    	{
    	    for(int s=0; s<OP_set_index; s++) { //for each set
    	    	op_set set = OP_set_list[s];
    	    	part P=OP_part_list[set->index];
    	    	if(P->is_partitioned != 1)
    	    	{
    	    	    printf("Unable to find mapping between primary set and %s \n",
    	    	    	P->set->name);
    	    	}
    	    }
    	    printf("Partitioning aborted !\n");
    	    exit(1);
    	}
    }

/*--STEP 2 - Create Imp/Export Lists for Migrating elements to new partitions-*/
    
    set_halo_list pe_list[OP_set_index]; //export list for each set 
    set_halo_list pi_list[OP_set_index]; //import list for each set
    
    //create partition export lists
    int* temp_list; int count, cap;
  
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=OP_set_list[s];
    	part p= OP_part_list[set->index];
      
    	//create a temporaty scratch space to hold export list for this set's 
    	//partition information
    	count = 0;cap = 1000;
    	temp_list = malloc(cap*sizeof(int));
    	if(temp_list == NULL) {
    	    printf(" partition_all -- error allocating memory: int* temp_list\n");
    	    exit(-1);
    	}
    
    	for(int i = 0; i < set->size; i++)
    	{
    	    if(p->elem_part[i] != my_rank)
    	    {
    	    	if(count>=cap)
    	    	{
    	    	    cap = cap*2;
    	    	    temp_list = realloc(temp_list, cap*sizeof(int));
    	    	    if(temp_list == NULL) {
    	    	    	printf(" partition_all -- error reallocating memory: int* temp_list\n");
    	    	    	exit(-1);
    	    	    }
    	    	}
    	    	temp_list[count++] = p->elem_part[i];
    	    	temp_list[count++] = i;//part.g_index[i];      	      
    	    }
    	}
    	//create partition info export list
    	pe_list[set->index] = (set_halo_list) malloc(sizeof(set_halo_list_core));
    	create_exp_list(set, temp_list, pe_list[set->index], count, comm_size, my_rank);
    	free(temp_list);      
    }
    
    
    //create partition import lists
   
    int *neighbors, *sizes;
    int ranks_size;
  
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=OP_set_list[s];
      
    	set_halo_list exp = pe_list[set->index];
      
    	//-----Discover neighbors-----
    	ranks_size = 0;
    	neighbors = malloc(comm_size*sizeof(int));
    	if(neighbors == NULL) {
    	    printf(" partition_all -- error allocating memory: int* neighbors\n");
    	    exit(-1);
    	}
    	sizes = malloc(comm_size*sizeof(int));
    	if(sizes == NULL) {
    	    printf(" partition_all -- error allocating memory: int* sizes\n");
    	    exit(-1);
    	}
      
	find_neighbors_set(exp, neighbors, sizes, &ranks_size, 
	    my_rank, comm_size, OP_PART_WORLD);	
    	MPI_Request request_send[exp->ranks_size];
      
    	int* rbuf; 
    	cap = 0; count = 0;
      
    	for(int i=0; i<exp->ranks_size; i++) {
    	    //printf("export from %d to %d set %10s, list of size %d \n",
    	    //my_rank,exp->ranks[i],set->name,exp->sizes[i]);
    	    int* sbuf = &exp->list[exp->disps[i]];
    	    MPI_Isend( sbuf,  exp->sizes[i],  MPI_INT, exp->ranks[i], 1,
    	    	OP_PART_WORLD, &request_send[i] );
    	}
      
    	for(int i=0; i< ranks_size; i++) cap = cap + sizes[i];
    	temp_list = malloc(cap*sizeof(int));
    	if(temp_list == NULL) {
    	    printf(" partition_all -- error allocating memory: int* temp_list\n");
    	    exit(-1);
    	}
      
    	for(int i=0; i<ranks_size; i++) {
    	    //printf("import from %d to %d set %10s, list of size %d\n",
    	    //neighbors[i], my_rank, set->name, sizes[i]);
    	    rbuf = malloc(sizes[i]*sizeof(int));
    	    
    	    MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], 1, OP_PART_WORLD,
    	    	MPI_STATUSES_IGNORE );
    	    memcpy(&temp_list[count],(void *)&rbuf[0],sizes[i]*sizeof(int));
    	    count = count + sizes[i];
    	    free(rbuf);
    	}
      
    	MPI_Waitall(exp->ranks_size,request_send, MPI_STATUSES_IGNORE );
    	pi_list[set->index] = (set_halo_list) malloc(sizeof(set_halo_list_core));
    	create_imp_list(set, temp_list, pi_list[set->index], count,
    	    neighbors, sizes, ranks_size, comm_size, my_rank);
    }
    

/*--STEP 3 - Perform Partitioning Data migration -----------------------------*/

    //data migration first ......  
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=OP_set_list[s];
    	
    	set_halo_list imp = pi_list[set->index];
    	set_halo_list exp = pe_list[set->index];
    	
    	MPI_Request request_send[exp->ranks_size];
    	    	
    	//migrate data defined on this set
    	for(int d=0; d<OP_dat_index; d++) { //for data array
    	    op_dat dat=OP_dat_list[d];
    	    
    	    if(compare_sets(dat->set,set)==1) //this data array is defines on this set
    	    {
    	    	
    	    	//prepare bits of the data array to be exported
    	    	char** sbuf = malloc(exp->ranks_size*sizeof(char *));
    	    	if(sbuf == NULL) {
    	    	    printf(" partition_all -- error allocating memory: char** sbuf\n");
    	    	    exit(-1);
    	    	}
    	
    	    	for(int i=0; i < exp->ranks_size; i++) {
    	    	    sbuf[i] = malloc(exp->sizes[i]*dat->size);
    	    	    if(sbuf[i] == NULL) {
    	    	    	printf(" partition_all -- error allocating memory: sbuf[i]\n");
    	    	    	exit(-1);
    	    	    }
    	    	
    	    	    for(int j = 0; j<exp->sizes[i]; j++)
    	    	    {
    	    	    	int index = exp->list[exp->disps[i]+j];
    	    	    	memcpy(&sbuf[i][j*dat->size],
    	    	    	    (void *)&dat->data[dat->size*(index)],dat->size);
    	    	    }
    	    	    //printf("export from %d to %d data %10s, number of elements of size %d | sending:\n ",
    	    	    //    my_rank,exp->ranks[i],dat->name,exp->sizes[i]);
    	    	    MPI_Isend(sbuf[i], dat->size*exp->sizes[i], MPI_CHAR, exp->ranks[i],
    	    	    	d, OP_PART_WORLD, &request_send[i]);
    	    	}
      	      
    	    	char *rbuf = malloc(dat->size*imp->size);
    	    	if(rbuf == NULL) {
    	    	    	printf(" partition_all -- error allocating memory: rbuf\n");
    	    	    	exit(-1);
    	    	}
    	    	
    	    	for(int i=0; i<imp->ranks_size; i++) {
    	    	    //printf("imported on to %d data %10s, number of elements of size %d | recieving:\n ",
    	    	    //	  my_rank, dat->name, imp->size);
    	    	    MPI_Recv(&rbuf[imp->disps[i]*dat->size],dat->size*imp->sizes[i], 
    	    	    	MPI_CHAR, imp->ranks[i], d, OP_PART_WORLD, MPI_STATUSES_IGNORE);
    	    	}
      	      
    	    	MPI_Waitall(exp->ranks_size,request_send, MPI_STATUSES_IGNORE );
    	    	for(int i=0; i < exp->ranks_size; i++) free(sbuf[i]); free(sbuf);
      	      
    	    	//delete the data entirs that has been sent and create a 
    	    	//modified data array
    	    	char* new_dat = malloc(dat->size*(set->size+imp->size));
    	    	if(new_dat == NULL) {
    	    	    	printf(" partition_all -- error allocating memory: char* new_dat\n");
    	    	    	exit(-1);
    	    	}
    	    	
    	    	int count = 0;
    	    	for(int i = 0; i < dat->set->size;i++)//iterate over old set size
    	    	{
    	    	    if(OP_part_list[set->index]->elem_part[i] == my_rank)
    	    	    {
    	    	    	memcpy(&new_dat[count*dat->size],(void *)&OP_dat_list[dat->index]->
    	    	    	    data[dat->size*i],dat->size);
    	    	    	count++;
    	    	    }      	      	            	      		  
    	    	}

    	    	memcpy(&new_dat[count*dat->size],(void *)rbuf,dat->size*imp->size);
    	    	count = count+imp->size;
    	    	new_dat = realloc(new_dat,dat->size*count);
    	    	if(new_dat == NULL & count != 0) {
    	    	    	printf(" partition_all -- error reallocating memory: char* new_dat\n");
    	    	    	exit(-1);
    	    	}
    	    	
    	    	free(rbuf);
      	           	      
    	    	free(OP_dat_list[dat->index]->data);   	  
    	    	OP_dat_list[dat->index]->data = new_dat;
    	    }  
    	}
    }
    
    //mapping tables second ......
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=OP_set_list[s];
      
    	set_halo_list imp = pi_list[set->index];
    	set_halo_list exp = pe_list[set->index];
      
    	MPI_Request request_send[exp->ranks_size];
           
    	//migrate mapping tables from this set
    	for(int m=0; m<OP_map_index; m++) { //for each maping table
    	    op_map map=OP_map_list[m];
      	  
    	    if(compare_sets(map->from,set)==1) //need to select mappings FROM this set
    	    {
    	    	//prepare bits of the mapping tables to be exported
    	    	int** sbuf = (int **)malloc(exp->ranks_size*sizeof(int *));
    	    	if(sbuf == NULL) {
    	    	    printf(" partition_all -- error allocating memory: int** sbuf\n");
    	    	    exit(-1);
    	    	}
    	    	//send mapping table entirs to relevant mpi processes
    	    	for(int i=0; i<exp->ranks_size; i++) {
    	    	    sbuf[i] = malloc(exp->sizes[i]*map->dim*sizeof(int));
    	    	    if(sbuf[i] == NULL) {
    	    	    	printf(" partition_all -- error allocating memory: sbuf[i]\n");
    	    	    	exit(-1);
    	    	    }
    	    	    for(int j = 0; j<exp->sizes[i]; j++)
    	    	    {
    	    	    	for(int p = 0; p< map->dim; p++)
    	    	    	{
    	    	    	    sbuf[i][j*map->dim+p] =
    	    	    	    map->map[map->dim*(exp->list[exp->disps[i]+j])+p];
    	    	    	}
    	    	    }      	      	  
    	    	    //printf("\n export from %d to %d map %10s, number of elements of size %d | sending:\n ",
    	    	    //    my_rank,exp->ranks[i],map->name,exp->sizes[i]);
    	    	    MPI_Isend(sbuf[i],  map->dim*exp->sizes[i],  MPI_INT, exp->ranks[i],
    	    	    	m, OP_PART_WORLD, &request_send[i]);
    	    	}
      	      
    	    	int *rbuf = malloc(map->dim*sizeof(int)*imp->size);
    	    	if(rbuf == NULL) {
    	    	    printf(" partition_all -- error allocating memory: int *rbuf \n");
    	    	    exit(-1);
    	    	}
      	      
    	    	//receive mapping table entirs from relevant mpi processes
    	    	for(int i=0; i < imp->ranks_size; i++) {
    	    	    //printf("\n imported on to %d map %10s, number of elements of size %d | recieving: ",
    	    	    //	  my_rank, map->name, imp->size);
    	    	    MPI_Recv(&rbuf[imp->disps[i]*map->dim], map->dim*imp->sizes[i], 
    	    	    	MPI_INT, imp->ranks[i], m, OP_PART_WORLD, MPI_STATUSES_IGNORE);
    	    	}
      	      
    	    	MPI_Waitall(exp->ranks_size,request_send, MPI_STATUSES_IGNORE );
    	    	for(int i=0; i < exp->ranks_size; i++) free(sbuf[i]); free(sbuf);
      	      
    	    	//delete the mapping table entirs that has been sent and create a 
    	    	//modified mapping table
    	    	int* new_map = malloc(sizeof(int)*(set->size+imp->size)*map->dim);
    	    	if(new_map == NULL) {
    	    	    printf(" partition_all -- error allocating memory: int* new_map \n");
    	    	    exit(-1);
    	    	}
    	    	
    	    	int count = 0;
    	    	for(int i = 0; i < map->from->size;i++)//iterate over old size of the maping table
    	    	{
    	    	    if(OP_part_list[map->from->index]->elem_part[i] == my_rank)
    	    	    {
    	    	    	memcpy(&new_map[count*map->dim],(void *)&OP_map_list[map->index]->
    	    	    	    map[map->dim*i],map->dim*sizeof(int));
    	    	    	count++;
    	    	    }
    	    	}
    	    	memcpy(&new_map[count*map->dim],(void *)rbuf,map->dim*sizeof(int)*imp->size);
    	    	count = count+imp->size;
    	    	new_map = realloc(new_map,sizeof(int)*count*map->dim);
    	    	if(new_map == NULL && count !=0) {
    	    	    printf(" partition_all -- error reallocating memory: int* new_map \n");
    	    	    exit(-1);
    	    	}
    	    	
    	    	free(rbuf);
    	    	free(OP_map_list[map->index]->map);
    	    	OP_map_list[map->index]->map = new_map;
    	    }   
    	}
    }

/*--STEP 4 - Update Partitioning Information and Sort Set Elements------------*/

    //update partition information and new local set sizes
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=OP_set_list[s];
      
    	part p = OP_part_list[set->index];
    	set_halo_list imp = pi_list[set->index];
      
    	int* new_g_index = malloc(sizeof(int)*(set->size+imp->size));
    	if(new_g_index == NULL ) {
    	    printf(" partition_all -- error allocating memory: int* new_g_index\n");
    	    exit(-1);
    	}
      
    	int count = 0;
    	for(int i = 0; i < set->size;i++)
    	{
    	    if(p->elem_part[i] == my_rank)
    	    {
    	    	new_g_index[count] = OP_part_list[set->index]->g_index[i];
    	    	count++;
    	    }      	  
    	}
      
    	for(int i = 0; i<imp->ranks_size; i++)
    	{
    	    for(int j = 0; j<imp->sizes[i]; j++)
    	    {
    	    	new_g_index[count] = get_global_index(imp->list[imp->disps[i]+j],
    	    	    imp->ranks[i], part_range[set->index], comm_size);
    	    	count++;
    	    }
    	}
      
    	new_g_index = realloc(new_g_index,sizeof(int)*count);
    	if(new_g_index == NULL && count != 0) {
    	    printf(" partition_all -- error reallocating memory: int* new_g_index\n");
    	    exit(-1);
    	}
    	
    	int* new_part = malloc(sizeof(int)*count);
    	if(new_part == NULL) {
    	    printf(" partition_all -- error allocating memory: int* new_part\n");
    	    exit(-1);
    	}
    	for(int i = 0; i< count; i++)new_part[i] = my_rank;
      
    	free(OP_part_list[set->index]->elem_part);
    	free(OP_part_list[set->index]->g_index);
      
    	OP_part_list[set->index]->elem_part = new_part;
    	OP_part_list[set->index]->g_index = new_g_index;   
      
    	OP_set_list[set->index]->size = count;   
    	OP_part_list[set->index]->set= OP_set_list[set->index];    

    }

   
  
    //re-set values in mapping tables
    for(int m=0; m<OP_map_index; m++) { //for each maping table
   	   	  op_map map=OP_map_list[m];
      	  
   	   	  OP_map_list[map->index]->from = OP_set_list[map->from->index];
   	   	  OP_map_list[map->index]->to = OP_set_list[map->to->index];
    }
  
    //re-set values in data arrays
    for(int d=0; d<OP_dat_index; d++) { //for data array
  	    	  op_dat dat=OP_dat_list[d];      	  
  	    	  OP_dat_list[dat->index]->set = OP_set_list[dat->set->index];      	  
    }

    //finally .... need to sort for each set, data on the set and mapping tables 
    //from this set accordiing to the OP_part_list[set.index]->g_index array values.  
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=OP_set_list[s];
                  
    	//first ... data on this set
    	for(int d=0; d<OP_dat_index; d++) { //for data array
    	    op_dat dat=OP_dat_list[d];
      	  
    	    if(compare_sets(dat->set,set) == 1)
    	    {
    	    	if(set->size > 0)
    	    	{
    	    	    int* temp = malloc(sizeof(int)*set->size);
    	    	    if(temp == NULL ) {
    	    	    	printf(" partition_all -- error allocating memory: int* temp\n");
    	    	    	exit(-1);
    	    	    }
    	    	    memcpy(temp, (void *)OP_part_list[set->index]->g_index, 
    	    	    	sizeof(int)*set->size);
    	    	    quickSort_dat(temp,OP_dat_list[dat->index]->data, 0, 
    	    	    	set->size-1, dat->size);
    	    	    free(temp);
    	    	}
    	    }
    	}
      
    	//second ... mapping tables
    	for(int m=0; m<OP_map_index; m++) { //for each maping table
    	    op_map map=OP_map_list[m];
    	    
    	    if(compare_sets(map->from,set) == 1)
    	    {
    	    	if(set->size > 0)
    	    	{
    	    	    int* temp = malloc(sizeof(int)*set->size);
    	    	    if(temp == NULL ) {
    	    	    	printf(" partition_all -- error allocating memory: int* temp\n");
    	    	    	exit(-1);
    	    	    }
    	    	    memcpy(temp, (void *)OP_part_list[set->index]->g_index, 
    	    	    	sizeof(int)*set->size);
    	    	    quickSort_map(temp,OP_map_list[map->index]->map, 0, 
    	    	    	set->size-1, map->dim); 
    	    	    free(temp);
    	    	}
    	    }      	  
    	}
    	if(set->size > 0) 
    	    quickSort(OP_part_list[set->index]->g_index, 0, set->size-1);
    }
  

/*--STEP 5 - Save old Partitioning Information -------------------------------*/
	     
    //allocate memory for old list
    OP_part_list_original = (part *)malloc(OP_set_index*sizeof(part));
    for(int s=0; s<OP_set_index; s++) { //for each set
    	op_set set=OP_set_list[s];
    	
    	part p = (part) malloc(sizeof(part_core));
    	p->g_index = OP_part_list[set->index]->g_index;
    	    	
    	OP_part_list_original[set->index] = p;    	
    }

    
/*--STEP 6 - Renumber mapping table entries with new partition's indexes------*/

    //update partition rage information 
    for(int i = 0; i<OP_set_index; i++)free(part_range[i]);
    get_part_range(part_range,my_rank,comm_size, OP_PART_WORLD);

    //find elements of the "to" set thats not in this local process
    for(int m=0; m<OP_map_index; m++) { //for each maping table
    	op_map map=OP_map_list[m];
      	  
      	int cap = 1000; int count = 0;
      	int* req_list = malloc(cap*sizeof(int));
      	if(req_list == NULL ) {
      	    printf(" partition_all -- error allocating memory: int* req_list\n");
      	    exit(-1);
      	}
    	    	    
      	for(int i = 0; i< map->from->size; i++)
      	{
      	    int local_index;//, global_index;
      	    for(int j=0; j<map->dim; j++)
      	    {
      	    	local_index = binary_search(OP_part_list[map->to->index]->g_index,
      	    	    map->map[i*map->dim+j], 0, map->to->size-1);
      	    	
      	    	if(count>=cap)
      	    	{
      	      	    cap = cap*2;
      	      	    req_list = realloc(req_list, cap*sizeof(int));
      	      	    if(req_list == NULL ) {
      	      	    	printf(" partition_all -- error reallocating memory: int* req_list\n");
      	      	    	exit(-1);
      	      	    }
      	      	}
      	      	  
      	      	if(local_index < 0) // not in this partition
      	      	{
      	      	    //store the global index of the element
      	      	    req_list[count++] = map->map[i*map->dim+j]; 
      	      	}
      	    }
      	}
      	//sort and remove duplicates
      	if(count > 0)
      	{
      	    //***point to check**//
      	    quickSort(req_list, 0, count-1);
      	    count = removeDups(req_list, count);
      	    req_list = realloc(req_list, count*sizeof(int));
      	    if(req_list == NULL && count != 0) {
      	    	printf(" partition_all -- error reallocating memory: int* req_list\n");
      	    	exit(-1);
      	    }
      	}
      	  
      	//do an allgather to findout how many elements that each process will 
      	//be asking partition information about
      	int recv_count[comm_size];
      	MPI_Allgather(&count, 1, MPI_INT, recv_count, 1, MPI_INT, OP_PART_WORLD);
      	  
      	//discover global size of these required elements
      	int g_count = 0;
      	for(int i = 0; i< comm_size; i++)g_count += recv_count[i];
      	  
      	//prepare for an allgatherv
      	int disp = 0;
      	int* displs = (int *)malloc(comm_size*sizeof(int));
      	if(displs == NULL) {
      	    printf(" partition_all -- error allocating memory: int* displs\n");
      	    exit(-1);
      	}
      	    
      	for(int i = 0; i<comm_size; i++)
      	{
      	    displs[i] =   disp;
      	    disp = disp + recv_count[i];
      	}
      	  
      	//allocate memory to hold the global indexes of elements requiring partition details
      	int *g_index = (int *)malloc(sizeof(int)*g_count);
      	if(g_index == NULL) {
      	    printf(" partition_all -- error allocating memory: int *g_index\n");
      	    exit(-1);
      	}
      	  
      	MPI_Allgatherv(req_list,count,MPI_INT, g_index,recv_count,displs,
      	    MPI_INT, OP_PART_WORLD);
      	free(req_list);
      	  
      	if(g_count > 0)
      	{
      	    quickSort(g_index, 0, g_count-1);
      	    g_count = removeDups(g_index, g_count);
      	    g_index = realloc(g_index, g_count*sizeof(int));
      	}
      	  
      	//printf("on rank %d map %s needs set %s : before g_count = %d\n", 
      	//    my_rank, map->name, map->to->name, g_count);
      	  
      	//go through the recieved global g_index array and see if any local element's
      	//partition details are requested by some foreign process
      	int *exp_index = (int *)malloc(sizeof(int)*g_count);
      	if(exp_index == NULL) {
      	    printf(" partition_all -- error allocating memory: int *exp_index\n");
      	    exit(-1);
      	}
      	
      	int *exp_g_index = (int *)malloc(sizeof(int)*g_count);
      	if(exp_g_index == NULL) {
      	    printf(" partition_all -- error allocating memory: int *exp_g_index\n");
      	    exit(-1);
      	}
      	
      	int exp_count = 0;
      	for(int i = 0; i<g_count; i++)
      	{
      	    int local_index = binary_search(OP_part_list[map->to->index]->g_index,
      	    	g_index[i],0,map->to->size-1);
      	    int global_index;
      	    if(local_index >= 0) 
      	    {
      	    	exp_g_index[exp_count] = g_index[i];
      	      	
      	    	global_index = get_global_index(local_index, my_rank,
      	    	    part_range[map->to->index], comm_size);
      	    	exp_index[exp_count++] = global_index;
      	    }
      	}
      	free(g_index);
      	
      	//realloc exp_index, exp_g_index
      	exp_index = (int *)realloc(exp_index,sizeof(int)*exp_count);
      	if(exp_index == NULL && exp_count != 0) {
      	    printf(" partition_all -- error reallocating memory: int *exp_index\n");
      	    exit(-1);
      	}
      	exp_g_index = (int *)realloc(exp_g_index,sizeof(int)*exp_count);
      	if(exp_g_index == NULL && exp_count != 0) {
      	    printf(" partition_all -- error reallocating memory: int *exp_g_index\n");
      	    exit(-1);
      	}
      	
      	//now export to every one these partition info with an all-to-all
      	MPI_Allgather(&exp_count, 1, MPI_INT, recv_count, 1, MPI_INT, OP_PART_WORLD);
      	disp = 0; free(displs);
      	displs = (int *)malloc(comm_size*sizeof(int));
      	if(displs == NULL ) {
      	    printf(" partition_all -- error allocating memory: int *displs\n");
      	    exit(-1);
      	}
      	
      	for(int i = 0; i<comm_size; i++)
      	{
      	    displs[i] =   disp;
      	    disp = disp + recv_count[i];
      	}
      	  
      	//allocate memory to hold the incomming partition details and allgatherv  
      	g_count = 0;
      	for(int i = 0; i< comm_size; i++)g_count += recv_count[i];
      	int *all_imp_index = (int *)malloc(sizeof(int)*g_count);
      	if(all_imp_index == NULL ) {
      	    printf(" partition_all -- error allocating memory: int *all_imp_index\n");
      	    exit(-1);
      	}
      	
      	g_index = malloc(sizeof(int)*g_count);
      	if(g_index == NULL ) {
      	    printf(" partition_all -- error allocating memory: int *g_index\n");
      	    exit(-1);
      	}
      	  
      	//printf("on rank %d map %s need set %s: After g_count = %d\n", 
      	//    my_rank, map.name,map.to.name,g_count);
      	  
      	MPI_Allgatherv(exp_g_index,exp_count,MPI_INT, g_index,recv_count,displs, 
      	    MPI_INT, OP_PART_WORLD);
      	
      	MPI_Allgatherv(exp_index,exp_count,MPI_INT, all_imp_index,recv_count,
      	    displs, MPI_INT, OP_PART_WORLD);
      	  
      	//sort all_imp_index according to g_index array
      	if(g_count > 0)quickSort_2(g_index, all_imp_index, 0, g_count-1); 
      	  
      	//now we hopefully we have all the informattion required to renumber this map
      	//so now, again go through each entry of this mapping table and renumber
      	for(int i = 0; i< map->from->size; i++)
      	{
      	    int local_index, global_index;
      	    for(int j=0; j < map->dim; j++)
      	    {
      	    	local_index = binary_search(OP_part_list[map->to->index]->g_index,
      	    	    map->map[i*map->dim+j], 0, map->to->size-1);
      	    	
      	    	if(local_index < 0) // not in this partition
      	      	{
      	      	    //need to search through g_index array
      	      	    int found = binary_search(g_index,map->map[i*map->dim+j],
      	      	    	0, g_count-1);
      	      	    if(found < 0) printf("Problem in renumbering\n");
      	      	    else
      	      	    {
      	      	    	OP_map_list[map->index]-> map[i*map->dim+j] = 
      	      	    	all_imp_index[found];
      	      	    }
      	      	}
      	      	else //in this partition
      	      	{
      	      	    global_index = get_global_index(local_index, my_rank,
      	      	    	part_range[map->to->index], comm_size);
      	      	    OP_map_list[map->index]->map[i*map->dim+j] = global_index;
      	      	}
      	    }
      	}
      	free(exp_index);free(exp_g_index);
      	free(g_index);
      	free(displs);
      	free(all_imp_index);
    } 
    
    
    //cleanup
    //destroy OP_part_list[], pe_list, pi_list    
    for(int s=0; s<OP_set_index; s++) { //for each set
      op_set set=OP_set_list[s];
      OP_part_list[set->index]->g_index = NULL;
      free(OP_part_list[set->index]->elem_part);
      
      free(pe_list[set->index]->ranks);free(pe_list[set->index]->disps);
      free(pe_list[set->index]->sizes);free(pe_list[set->index]->list);
      
      free(pi_list[set->index]->ranks);free(pi_list[set->index]->disps);
      free(pi_list[set->index]->sizes);free(pi_list[set->index]->list);
    }
    for(int i = 0; i<OP_set_index; i++)free(part_range[i]);free(part_range);
}

    


/**------------Partition A Primary Set Using XYZ Geometry Data---------------**/
void op_partition_geom(op_dat coords, int g_nnode)
{
    //declare timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    double time;
    double max_time;
    
    op_timers(&cpu_t1, &wall_t1); //timer start for partitioning
  
    //create new communicator for partitioning
    int my_rank, comm_size;
    MPI_Comm_dup(MPI_COMM_WORLD, &OP_PART_WORLD);
    MPI_Comm_rank(OP_PART_WORLD, &my_rank);
    MPI_Comm_size(OP_PART_WORLD, &comm_size);

/*--STEP 0 - initialise partitioning data stauctures with the current block 
    partitioning information */
    
    // Compute global partition range information for each set
    int** part_range = (int **)malloc(OP_set_index*sizeof(int*));
    if(part_range == NULL) {
    	printf(" op_partition_geom -- error allocating memory: int** part_range\n");
    	exit(-1);
    }
    get_part_range(part_range,my_rank,comm_size, OP_PART_WORLD);
    
    //allocate memory for list
    OP_part_list = (part *)malloc(OP_set_index*sizeof(part));
    
    for(int s=0; s<OP_set_index; s++) { //for each set
      op_set set=OP_set_list[s];
      //printf("set %s size = %d\n", set.name, set.size);
      int *g_index = (int *)malloc(sizeof(int)*set->size);
      if(g_index == NULL) {
      	  printf(" op_partition_geom -- error allocating memory: int *g_index\n");
      	  exit(-1);
      }
      
      for(int i = 0; i< set->size; i++)
      	  g_index[i] = get_global_index(i,my_rank, part_range[set->index],comm_size);      
      decl_partition(set, g_index, NULL); 
    }
  
/*--STEP 1 - Partition Nodes (primary set) using Coordinates (1D,2D or 3D)----*/

    //Setup data structures for ParMetis PartGeom
    idxtype *vtxdist = (idxtype *)malloc(sizeof(idxtype)*(comm_size+1));
    if(vtxdist == NULL) {
      	  printf(" op_partition_geom -- error allocating memory: idxtype *vtxdist\n");
      	  exit(-1);
    }
      
    idxtype *partition = (idxtype *)malloc(sizeof(idxtype)*coords->set->size);
    if(partition == NULL) {
      	  printf(" op_partition_geom -- error allocating memory: idxtype *partition\n");
      	  exit(-1);
    }
    
    int ndims = coords->dim;
    float* xyz;
    
    // Create ParMetis compatible coordinates array 
    //- i.e. coordinates should be floats
    if(ndims == 3 || ndims == 2 || ndims == 1)
    {
    	xyz = (float* )malloc(coords->set->size*coords->dim*sizeof(float));
    	if(xyz == NULL) {
    	    printf(" op_partition_geom -- error allocating memory: float* xyz\n");
    	    exit(-1);
      	}
    	size_t mult = coords->size/coords->dim;
    	for(int i = 0;i < coords->set->size;i++)
    	{
    	    double temp;
    	    for(int e = 0; e < coords->dim;e++)
    	    {
    	    	memcpy(&temp, (void *)&(OP_dat_list[coords->index]->
    	    	    data[(i*coords->dim+e)*mult]), mult);
    	    	xyz[i*coords->dim + e] = (float)temp;
    	    }
    	}
    }
    else
    {
    	printf("Dimensions of Coordinate array not one of 3D,2D or 1D\n");
    	printf("Not supported by ParMetis - Indicate correct coordinates array\n");
    	exit(1);
    }
  
    for(int i=0; i<comm_size; i++)
    {
    	vtxdist[i] = part_range[0][2*i];//nodes have index 0
    }
    vtxdist[comm_size] = g_nnode;
  
    //use xyz coordinates to feed into ParMETIS_V3_PartGeom
    ParMETIS_V3_PartGeom(vtxdist, &ndims, xyz, partition, &OP_PART_WORLD);
    free(xyz);free(vtxdist);
  
    //initialise primary set as partitioned
    OP_part_list[coords->set->index]->elem_part= partition;
    OP_part_list[coords->set->index]->is_partitioned = 1;
  
    
    //partition all other sets, migrate data and renumber mapping tables
    partition_all(coords->set, my_rank, comm_size);
    
    //cleanup
    for(int i = 0; i<OP_set_index; i++)free(part_range[i]);free(part_range);
    
    op_timers(&cpu_t2, &wall_t2);  //timer stop for partitioning
    //printf time for partitioning
    time = wall_t2-wall_t1;
    MPI_Reduce(&time,&max_time,1,MPI_DOUBLE, MPI_MAX,0, OP_PART_WORLD);
    MPI_Comm_free(&OP_PART_WORLD);  
    if(my_rank==0)printf("Max total geometric partitioning time = %lf\n",max_time);    
}


//destroy OP_part_list[]    
void op_partition_destroy()
{
    for(int s=0; s<OP_set_index; s++) { //for each set
      op_set set=OP_set_list[s];
      free(OP_part_list_original[set->index]->g_index);
      //free(OP_part_list_original[set->index]->elem_part);
    }     
}
