
//
// header file (includes op_lib_core.h and various system header files)
//

#include "op_lib.h"

// routines called by user code and kernels

void op_init(int argc, char **argv, int diags){
  op_init_core(argc, argv, diags);
}

void op_fetch_data(op_dat data) {}
void op_decl_const_char(int, char const*, int, char*, char const*){}

void op_exit(){
   for(int ip=0; ip<OP_plan_index; ip++) {
    for (int m=0; m<OP_plans[ip].nargs; m++)
      if (OP_plans[ip].maps[m] != NULL)
        free(OP_plans[ip].maps[m]);
    for (int m=0; m<OP_plans[ip].ninds; m++)
      free(OP_plans[ip].ind_maps[m]);
    free(OP_plans[ip].ind_offs);
    free(OP_plans[ip].ind_sizes);
    free(OP_plans[ip].nthrcol);
    free(OP_plans[ip].thrcol);
    free(OP_plans[ip].offset);
    free(OP_plans[ip].nelems);
    free(OP_plans[ip].blkmap);
  }

  op_exit_core();
}


// routines called by op_lib_core.cpp

extern "C"
void op_mvHostToDevice(void **map, int size) {}

extern "C"
void op_cpHostToDevice(void **dat_h, void **dat_d, int size) {}


