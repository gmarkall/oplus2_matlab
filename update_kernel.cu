// 
// auto-generated by op2.m on 07-Apr-2010 12:41:20 
//

// user function                                            
                                                            
__device__                                                  
#include "update.h"                                         
                                                            
                                                            
// CUDA kernel function                                     
                                                            
__global__ void op_cuda_update( float  *arg0,               
                                float  *arg1,               
                                float  *arg2,               
                                int set_size ) {            
                                                            
  // process set elements                                   
                                                            
  for (int n=threadIdx.x+blockIdx.x*blockDim.x;             
       n<set_size; n+=blockDim.x*gridDim.x) {               
                                                            
    // user-supplied kernel call                            
                                                            
    update( arg0+n*1,                                       
            arg1+n*1,                                       
            arg2+n*1 );                                     
  }                                                         
}                                                           
                                                            
                                                            
// host stub function                                       
                                                            
extern "C"                                                  
void op_par_loop_update(char const * name, op_set set,      
     op_dat arg0, int arg0idx, op_ptr arg0ptr, int arg0dim, 
          op_datatype arg0typ,           op_access arg0acc, 
     op_dat arg1, int arg1idx, op_ptr arg1ptr, int arg1dim, 
          op_datatype arg1typ,           op_access arg1acc, 
     op_dat arg2, int arg2idx, op_ptr arg2ptr, int arg2dim, 
          op_datatype arg2typ,           op_access arg2acc){
                                                            
  if (OP_DIAGS>1) {                                         
    printf(" kernel routine w/o indirection:  update \n");  
  }                                                         
                                                            
  // execute plan                                           
                                                            
  op_cuda_update<<<100,64>>>( (float  *) arg0.dat_d,        
                              (float  *) arg1.dat_d,        
                              (float  *) arg2.dat_d,        
                              set.size );                   
                                                            
  cutilCheckMsg("op_cuda_update execution failed\n");       
}                                                           
                                                            
