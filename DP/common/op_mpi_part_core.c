//parmetis header
#include <parmetis.h>

#include "op_mpi_part_core.h"

//
//MPI Communicator for partitioning
//
MPI_Comm OP_PART_WORLD;


/**----------------MPI partitioning related global variables ----------------**/
part* OP_part_list; 


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
    
    //create_list(map.to, temp_list, pe_list, count,
    	//neighbors, sizes, ranks_size, comm_size, my_rank);
    
    
    
    //OP_part_info_list[map->from->index]->part = partition;
    OP_part_list[map->from->index]->is_partitioned = 1;
    
    //cleanup
    free(pi_list->list);free(pi_list->ranks);free(pi_list->sizes);
    free(pi_list->disps);free(pi_list);
    
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
    	printf(" op_partition_geom -- error allocating memory: int** part_range\n");
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

/*--STEP 3 - Perform Partitioning Data migration -----------------------------*/

/*--STEP 4 - Update Partitioning Information and Sort Set Elements------------*/

/*--STEP 5 - Save old Partitioning Information -------------------------------*/
	     //create mapping from old to new

/*--STEP 6 - Renumber mapping table entries with new partition's indexes------*/

    //cleanup
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
      free(OP_part_list[set->index]->g_index);
      free(OP_part_list[set->index]->elem_part);
    }   
    
}
