// 
// auto-generated by op2.m on 04-Aug-2010 11:26:18 
//

// user function                                                                  
                                                                                  
__device__                                                                        
#include "update.h"                                                               
                                                                                  
                                                                                  
// CUDA kernel function                                                           
                                                                                  
__global__ void op_cuda_update(                                                   
  float *arg0,                                                                    
  float *arg1,                                                                    
  float *arg2,                                                                    
  float *arg3,                                                                    
  float *arg4,                                                                    
  int   set_size ) {                                                              
                                                                                  
  float arg4_l[1];                                                                
  for (int d=0; d<1; d++) arg4_l[d]=ZERO_float;                                   
                                                                                  
  // process set elements                                                         
                                                                                  
  for (int n=threadIdx.x+blockIdx.x*blockDim.x;                                   
       n<set_size; n+=blockDim.x*gridDim.x) {                                     
                                                                                  
      // user-supplied kernel call                                                
                                                                                  
      update( arg0+n*4,                                                           
              arg1+n*4,                                                           
              arg2+n*4,                                                           
              arg3+n*1,                                                           
              arg4_l );                                                           
  }                                                                               
                                                                                  
  // global reductions                                                            
                                                                                  
  for(int d=0; d<1; d++) op_reduction<OP_INC>(&arg4[d],arg4_l[d]);                
}                                                                                 
                                                                                  
                                                                                  
// host stub function                                                             
                                                                                  
void op_par_loop_update(char const *name, op_set set,                             
  op_dat arg0, int idx0, op_ptr ptr0, int dim0, char const *typ0, op_access acc0, 
  op_dat arg1, int idx1, op_ptr ptr1, int dim1, char const *typ1, op_access acc1, 
  op_dat arg2, int idx2, op_ptr ptr2, int dim2, char const *typ2, op_access acc2, 
  op_dat arg3, int idx3, op_ptr ptr3, int dim3, char const *typ3, op_access acc3, 
  float *arg4h,int idx4, op_ptr ptr4, int dim4, char const *typ4, op_access acc4){
                                                                                  
  op_dat arg4 = {{0,0,"null"},0,0,0,(char *)arg4h,NULL,"float","gbl"};            
                                                                                  
  if (OP_DIAGS>1) {                                                               
    printf(" kernel routine w/o indirection:  update \n");                        
  }                                                                               
                                                                                  
  // transfer global reduction data to GPU                                        
                                                                                  
  int reduct_bytes = 0;                                                           
  int reduct_size  = 0;                                                           
  reduct_bytes += ROUND_UP(1*sizeof(float));                                      
  reduct_size   = MAX(reduct_size,sizeof(float));                                 
                                                                                  
  reallocReductArrays(reduct_bytes);                                              
                                                                                  
  reduct_bytes = 0;                                                               
  arg4.dat   = OP_reduct_h + reduct_bytes;                                        
  arg4.dat_d = OP_reduct_d + reduct_bytes;                                        
  for (int d=0; d<1; d++) ((float *)arg4.dat)[d] = arg4h[d];                      
  reduct_bytes += ROUND_UP(1*sizeof(float));                                      
                                                                                  
  mvReductArraysToDevice(reduct_bytes);                                           
                                                                                  
  // execute plan                                                                 
                                                                                  
  int nshared = reduct_size*64/2;                                                 
                                                                                  
  op_cuda_update<<<100,64,nshared>>>( (float *) arg0.dat_d,                       
                                      (float *) arg1.dat_d,                       
                                      (float *) arg2.dat_d,                       
                                      (float *) arg3.dat_d,                       
                                      (float *) arg4.dat_d,                       
                                      set.size );                                 
                                                                                  
  cutilCheckMsg("op_cuda_update execution failed\n");                             
                                                                                  
  // transfer global reduction data back to CPU                                   
                                                                                  
  mvReductArraysToHost(reduct_bytes);                                             
                                                                                  
  for (int d=0; d<1; d++) arg4h[d] = ((float *)arg4.dat)[d];                      
}                                                                                 
                                                                                  
