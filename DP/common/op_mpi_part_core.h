

/**-------------------Data structures related to partitioning----------------**/

//struct to hold the partition information for each set
typedef struct
{
  op_set set; //set to which this partition information blongs to 
  int* g_index; //global index of each element held in this MPI process
  int* elem_part; //partition to which each element belongs
  int is_partitioned; //indicates if this set is partitioned 1 if partitioned 0 if not
} part_core;

typedef part_core* part;
