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
 
#include <stdlib.h>                                                                      
#include <stdio.h>                                                                       
#include <string.h>                                                                      
#include <math.h>
     
//                                                                                       
// templates                                                                             
//                                                                                       

template < class T, class D >                                                            
void arg_check(op_set set, int m, const char *name, D *arg, int idx, op_map map,         
               int dim, const char *typ, op_access acc, T **p_arg, int *ninds){          
  if (idx != -1) {                                                                       
    printf("error: arg %d in kernel \"%s\"\n",m,name);                                   
    printf("invalid index, should be -1 for constant arg\n");                            
    exit(1);                                                                             
  }                                                                                      
  if (type_error(arg,typ)) {                                                             
    printf("error: arg %d in kernel \"%s\"\n",m,name);                                   
    printf("data type does not match type of input dataset\n");                          
    exit(1);                                                                             
  }                                                                                      
  if (type_error(*p_arg,typ)) {                                                          
    printf("error: arg %d in kernel \"%s\"\n",m,name);                                   
    printf("data type does not match type of function argument\n");                      
    exit(1);                                                                             
  }                                                                                      
  if (dim <= 0) {                                                                        
    printf("error: arg %d in kernel \"%s\"\n",m,name);                                   
    printf("dimension must be strictly positive\n");                                     
    exit(1);                                                                             
  }                                                                                      
}                                                                                        
                                                                                         
template < class T >                                                                     
void arg_check(op_set set, int m, const char *name, op_dat arg, int idx, op_map map,     
               int dim, const char *typ, op_access acc, T **p_arg, int *ninds){          
  if (idx>=0) (*ninds)++;                                                                
                                                                                         
  if (map == NULL) {                                                                     
    if (idx != -1 && map->from->size != 0) {      //CHANGED                                                                      
      printf("error: arg %d in kernel \"%s\"\n",m,name);                                 
      printf("invalid index, should be -1 for identity mapping\n");                      
      exit(1);                                                                           
    }                                                                                    
  }                                                                                      
  else {                                                                                 
    if (set != map->from || arg->set != map->to) {                                       
      printf("error: arg %d in kernel \"%s\"\n",m,name);                                 
      printf("invalid pointer\n");                                                       
      exit(1);                                                                           
    }                                                                                    
    if (idx < 0 || idx >= map->dim) {                                                    
      printf("error: arg %d in kernel \"%s\"\n",m,name);                                 
      printf("invalid pointer index\n");                                                 
      exit(1);                                                                           
    }                                                                                    
  }                                                                                      
  if (strcmp(arg->type,typ)) {                                                           
    printf("error: arg %d in kernel \"%s\"\n",m,name);                                   
    printf("data type does not match type of input dataset\n");                          
    exit(1);                                                                             
  }                                                                                      
  if (type_error(*p_arg,typ)) {                                                          
    printf("error: arg %d in kernel \"%s\"\n",m,name);                                   
    printf("data type does not match type of function argument\n");                      
    exit(1);                                                                             
  }                                                                                      
  if (arg->dim != dim) {                                                                 
    printf("error: arg %d in kernel \"%s\"\n",m,name);                                   
    printf("dimension does not match input dataset\n");                                  
    exit(1);                                                                             
  }                                                                                      
}                                                                                        
                                                                                         
template < class T, class D >                                                            
void arg_set(int n,D *arg,int idx,op_map map,int dim,                                    
             const char *typ,op_access acc,T **p_arg){                                   
  *p_arg = (T *) arg;                                                                    
}                                                                                        
                                                                                         
template < class T >                                                                     
void arg_set(int n,op_dat arg,int idx,op_map map,int dim,                                
             const char *typ,op_access acc,T **p_arg){                                   
  int n2;                                                                                
  if (map==NULL)                    // identity mapping                                  
    n2 = n;                                                                              
  else                              // standard pointers                                 
    n2 = map->map[idx+n*map->dim];                                                       
                                                                                         
  *p_arg = (T *)(arg->data + n2*arg->size);                                               
}


//blank out arguments with global reductions for indirect loop halo executions
template <class D>
D* blank_arg(D* arg, op_map map, op_access acc)
{
  /**check if returning NULL is ok for this**/
  double junck = 0.0;
  if(/*map->index == -2 && */acc != OP_READ)//this argument is OP_GBL and OP_INC or OP_MAX/MIN
  {   
      return NULL;//(D *)&junck; 
  }
  else 
  {
      return arg;
  }
}


//                                                                                       
// op_par_loop routine for 2 arguments    - save_soln                                               
//                                                                                       
            
template < class T0, class T1,                                                           
           class D0, class D1 >                                                          
void op_par_loop_save_soln(void (*kernel)( T0*, T1* ),                                             
  char const * name, op_set set,                                                         
  D0  arg0 ,int idx0 ,op_map map0 ,int dim0 ,const char *typ0 ,op_access acc0 ,          
  D1  arg1 ,int idx1 ,op_map map1 ,int dim1 ,const char *typ1 ,op_access acc1 ) 
{
    char *p_arg0, *p_arg1;  
    int ninds=0;
    int exec_length = 0;
    
    // consistency checks
    
    if (OP_diags>0) {
    	arg_check(set,0 ,name,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 ,&ninds);
    	arg_check(set,1 ,name,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 ,&ninds);
    }
    
    if (OP_diags>2) {
    	if (ninds==0)
    	    printf(" kernel routine w/o indirection:  %s \n",name);
    	else
    	    printf(" kernel routine with indirection: %s \n",name);
    }
    
    // initialise timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    op_timers(&cpu_t1, &wall_t1); 
    
    if(idx0 != -1 || idx1 != -1)//indirect loop
    {
    	//for each indirect data set
    	exchange_halo(set, arg0, acc0, idx0);
    	exchange_halo(set, arg1, acc1, idx1);
    	
    	//for all indirect dataset access with OP_READ
    	if(acc0 == OP_READ && acc1 == OP_READ   ) exec_length = set->size;
    	else  exec_length = set->size + OP_import_sets_list[set->index]->size;
    }
    else //direct loop
    {
    	exec_length = set->size;
    }
    
    // loop over set elements
    //(1) over owned partition
    for (int n=0; n<set->size; n++) {
    	arg_set(n,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 );
    	arg_set(n,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 );
    	
    	// call kernel function, passing in pointers to data
    	kernel( (T0 *)p_arg0,  (T1 *)p_arg1 );
    }
    
    //(2) over exec halo (blank out global parameters to avoid double counting)
    for (int n=set->size; n<exec_length; n++) {
    	arg_set(n,*(blank_arg(&arg0,map0,acc0)),idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 );
    	arg_set(n,*(blank_arg(&arg1,map1,acc1)),idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 );
    	
    	// call kernel function, passing in pointers to exec halo data
    	kernel( (T0 *)p_arg0,  (T1 *)p_arg1 );
    }
  
    //set dirty bit on direct/indirect datasets with access OP_INC,OP_WRITE, OP_RW
    set_dirtybit(arg0, acc0);
    set_dirtybit(arg1, acc1);

  
    //performe any global operations
    //--none
    
  
    //update timer record
    op_timers(&cpu_t2, &wall_t2);
    
    op_timing_realloc(0);
    if(OP_kernels[0].count==0)
    	OP_kernels[0].name     = name;
    OP_kernels[0].count    += 1;
    OP_kernels[0].time     += wall_t2 - wall_t1;                               
} 



//                                                                                       
// op_par_loop routine for 5 arguments     - update                                               
//                                                                                      
template < class T0, class T1, class T2, class T3, class T4,                             
           class D0, class D1, class D2, class D3, class D4 >                            
void op_par_loop_update(void (*kernel)( T0*, T1*, T2*, T3*, T4* ),                            
  char const * name,  op_set set,                                                         
  D0  arg0 ,int idx0 ,op_map map0 ,int dim0 ,const char *typ0 ,op_access acc0 ,          
  D1  arg1 ,int idx1 ,op_map map1 ,int dim1 ,const char *typ1 ,op_access acc1 ,          
  D2  arg2 ,int idx2 ,op_map map2 ,int dim2 ,const char *typ2 ,op_access acc2 ,          
  D3  arg3 ,int idx3 ,op_map map3 ,int dim3 ,const char *typ3 ,op_access acc3 ,          
  D4  arg4 ,int idx4 ,op_map map4 ,int dim4 ,const char *typ4 ,op_access acc4 ){          
                                                                                         
  char *p_arg0, *p_arg1, *p_arg2, *p_arg3, *p_arg4;                                                                          
                                                                                         
  int ninds=0;
  int exec_length = 0;                                                                          
                     
  
  // consistency checks                                                                  
                                                                                         
  if (OP_diags>0) {                                                                      
    arg_check(set,0 ,name,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 ,&ninds);  
    arg_check(set,1 ,name,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 ,&ninds);  
    arg_check(set,2 ,name,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 ,&ninds);  
    arg_check(set,3 ,name,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 ,&ninds);  
    arg_check(set,4 ,name,arg4 ,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 ,&ninds);   
  }                                                                                      
                                                                                         
 if (OP_diags>2) {                                                                      
    if (ninds==0)                                                                        
      printf(" kernel routine w/o indirection:  %s \n",name);                            
    else                                                                                 
      printf(" kernel routine with indirection: %s \n",name);                            
  }
  
  // initialise timers                                                            
  double cpu_t1, cpu_t2, wall_t1, wall_t2;                                        
  op_timers(&cpu_t1, &wall_t1); 
 
  if(idx0 != -1 || idx1 != -1 || idx2 != -1 || 
     idx3 != -1 || idx4 != -1 )//indirect loop
  {
      //for each indirect data set
      exchange_halo(set, arg0, acc0, idx0); 
      exchange_halo(set, arg1, acc1, idx1);
      exchange_halo(set, arg2, acc2, idx2);
      exchange_halo(set, arg3, acc3, idx3);
      exchange_halo_void(set, arg4, acc4, idx4);

             
      //for all indirect dataset access with OP_READ
      if(acc0 == OP_READ && acc1 == OP_READ && acc2 == OP_READ && acc3 == OP_READ &&
      	 acc4 == OP_READ ) exec_length = set->size;
      else  exec_length = set->size + OP_import_sets_list[set->index]->size;
  }
  else //direct loop
  {
      exec_length = set->size;
  }
                                                                                         
  // loop over set elements 
  //(1) over owned partition 
  for (int n=0; n < set->size; n++) {                                                       
    arg_set(n,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 );                     
    arg_set(n,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 );                     
    arg_set(n,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 );                     
    arg_set(n,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 );                     
    arg_set(n,arg4 ,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 );     
    
    kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,  (T4 *)p_arg4);                                                              
  }             
  
  //***arg 4 id a OP_GBL reduction - need to stop double counting***
  double junck = 0.0;
  
  //(2) over exec halo (blank out global parameters to avoid double counting)
  for (int n=set->size; n<exec_length; n++) { 
    arg_set(n,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 );                     
    arg_set(n,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 );                     
    arg_set(n,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 );                     
    arg_set(n,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 );                     
    arg_set(n,&junck,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 );  
      
      // call kernel function, passing in pointers to exec halo data  
      kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,  (T4 *)p_arg4);    
  }
  
  //set dirty bit on direct/indirect datasets with access OP_INC,OP_WRITE, OP_RW
  set_dirtybit(arg0, acc0);
  set_dirtybit(arg1, acc1);
  set_dirtybit(arg2, acc2);
  set_dirtybit(arg3, acc3);
  set_dirtybit_void(arg4, acc4); //no dirty bit 
  
  //performe any global operations
  global_reduce_double(arg4, acc4);
    
  //update timer record
  op_timers(&cpu_t2, &wall_t2);
  op_timing_realloc(4);
  if(OP_kernels[4].count==0)
      OP_kernels[4].name     = name;
  OP_kernels[4].count    += 1;
  OP_kernels[4].time     += wall_t2 - wall_t1;  
}  



                                                                                         
//                                                                                       
// op_par_loop routine for 6 arguments       - adt                                             
//                                                                                       
                                                                                         
template < class T0, class T1, class T2, class T3, class T4,                             
           class T5,                                                                     
           class D0, class D1, class D2, class D3, class D4,                             
           class D5 >                                                                    
void op_par_loop_adt(void (*kernel)( T0*, T1*, T2*, T3*, T4*,                              
                                    T5* ),                                               
  char const * name, op_set set,                                                         
  D0  arg0 ,int idx0 ,op_map map0 ,int dim0 ,const char *typ0 ,op_access acc0 ,          
  D1  arg1 ,int idx1 ,op_map map1 ,int dim1 ,const char *typ1 ,op_access acc1 ,          
  D2  arg2 ,int idx2 ,op_map map2 ,int dim2 ,const char *typ2 ,op_access acc2 ,          
  D3  arg3 ,int idx3 ,op_map map3 ,int dim3 ,const char *typ3 ,op_access acc3 ,          
  D4  arg4 ,int idx4 ,op_map map4 ,int dim4 ,const char *typ4 ,op_access acc4 ,          
  D5  arg5 ,int idx5 ,op_map map5 ,int dim5 ,const char *typ5 ,op_access acc5 ){         
                                                                                         
  char *p_arg0, *p_arg1, *p_arg2, *p_arg3, *p_arg4,                                      
       *p_arg5;                                                                          
                                                                                         
  int ninds=0;
  int exec_length = 0;                                                                          
                     
  
  // consistency checks                                                                  
                                                                                         
  if (OP_diags>0) {                                                                      
    arg_check(set,0 ,name,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 ,&ninds);  
    arg_check(set,1 ,name,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 ,&ninds);  
    arg_check(set,2 ,name,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 ,&ninds);  
    arg_check(set,3 ,name,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 ,&ninds);  
    arg_check(set,4 ,name,arg4 ,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 ,&ninds);  
    arg_check(set,5 ,name,arg5 ,idx5 ,map5 ,dim5 ,typ5 ,acc5 ,(T5  **)&p_arg5 ,&ninds);  
  }                                                                                      
                                                                                         
 if (OP_diags>2) {                                                                      
    if (ninds==0)                                                                        
      printf(" kernel routine w/o indirection:  %s \n",name);                            
    else                                                                                 
      printf(" kernel routine with indirection: %s \n",name);                            
  }
  
  // initialise timers                                                            
  double cpu_t1, cpu_t2, wall_t1, wall_t2;                                        
  op_timers(&cpu_t1, &wall_t1); 
 
  if(idx0 != -1 || idx1 != -1 || idx2 != -1 || 
     idx3 != -1 || idx4 != -1 || idx5 != -1)//indirect loop
  {
      //for each indirect data set
      exchange_halo(set, arg0, acc0, idx0); 
      exchange_halo(set, arg1, acc1, idx1);
      exchange_halo(set, arg2, acc2, idx2);
      exchange_halo(set, arg3, acc3, idx3);
      exchange_halo(set, arg4, acc4, idx4);
      exchange_halo(set, arg5, acc5, idx5);
             
      //for all indirect dataset access with OP_READ
      if(acc0 == OP_READ && acc1 == OP_READ && acc2 == OP_READ && acc3 == OP_READ &&
      	 acc4 == OP_READ && acc5 == OP_READ ) exec_length = set->size;
      else  exec_length = set->size + OP_import_sets_list[set->index]->size;
  }
  else //direct loop
  {
      exec_length = set->size;
  }
                                                                                         
  // loop over set elements 
  //(1) over owned partition 
  for (int n=0; n < set->size; n++) {                                                       
    arg_set(n,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 );                     
    arg_set(n,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 );                     
    arg_set(n,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 );                     
    arg_set(n,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 );                     
    arg_set(n,arg4 ,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 );                     
    arg_set(n,arg5 ,idx5 ,map5 ,dim5 ,typ5 ,acc5 ,(T5  **)&p_arg5 );                 
    
    kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,  (T4 *)p_arg4,    
            (T5 *)p_arg5 );                                                              
  }             
  
  
  //(2) over exec halo (blank out global parameters to avoid double counting)
  for (int n=set->size; n<exec_length; n++) { 
    arg_set(n,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 );                     
    arg_set(n,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 );                     
    arg_set(n,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 );                     
    arg_set(n,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 );                     
    arg_set(n,arg4 ,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 );                     
    arg_set(n,arg5 ,idx5 ,map5 ,dim5 ,typ5 ,acc5 ,(T5  **)&p_arg5 );   
      
      // call kernel function, passing in pointers to exec halo data  
      kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,  (T4 *)p_arg4,    
            (T5 *)p_arg5 );    
  }
  
  //set dirty bit on direct/indirect datasets with access OP_INC,OP_WRITE, OP_RW
  set_dirtybit(arg0, acc0);
  set_dirtybit(arg1, acc1);
  set_dirtybit(arg2, acc2);
  set_dirtybit(arg3, acc3);
  set_dirtybit(arg4, acc4);
  set_dirtybit(arg5, acc5);

  
  //performe any global operations
  //--none

  //update timer record
  op_timers(&cpu_t2, &wall_t2);
  op_timing_realloc(1);
  if(OP_kernels[1].count==0)
      OP_kernels[1].name     = name;
  OP_kernels[1].count    += 1;
  OP_kernels[1].time     += wall_t2 - wall_t1; 
  
}                                                                                        
 

//                                                                                       
// op_par_loop routine for 6 arguments       - bres_calc                                             
//                                                                                       
template < class T0, class T1, class T2, class T3, class T4,                             
           class T5,                                                                     
           class D0, class D1, class D2, class D3, class D4,                             
           class D5 >                                                                    
void op_par_loop_bres_calc(void (*kernel)( T0*, T1*, T2*, T3*, T4*,                              
                                    T5* ),                                               
  char const * name, op_set set,                                                         
  D0  arg0 ,int idx0 ,op_map map0 ,int dim0 ,const char *typ0 ,op_access acc0 ,          
  D1  arg1 ,int idx1 ,op_map map1 ,int dim1 ,const char *typ1 ,op_access acc1 ,          
  D2  arg2 ,int idx2 ,op_map map2 ,int dim2 ,const char *typ2 ,op_access acc2 ,          
  D3  arg3 ,int idx3 ,op_map map3 ,int dim3 ,const char *typ3 ,op_access acc3 ,          
  D4  arg4 ,int idx4 ,op_map map4 ,int dim4 ,const char *typ4 ,op_access acc4 ,          
  D5  arg5 ,int idx5 ,op_map map5 ,int dim5 ,const char *typ5 ,op_access acc5 ){         
                                                                                         
  char *p_arg0, *p_arg1, *p_arg2, *p_arg3, *p_arg4,                                      
       *p_arg5;                                                                          
                                                                                         
  int ninds=0;
  int exec_length = 0;                                                                          
                     
  
  // consistency checks                                                                  
                                                                                         
  if (OP_diags>0) {                                                                      
    arg_check(set,0 ,name,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 ,&ninds);  
    arg_check(set,1 ,name,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 ,&ninds);  
    arg_check(set,2 ,name,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 ,&ninds);  
    arg_check(set,3 ,name,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 ,&ninds);  
    arg_check(set,4 ,name,arg4 ,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 ,&ninds);  
    arg_check(set,5 ,name,arg5 ,idx5 ,map5 ,dim5 ,typ5 ,acc5 ,(T5  **)&p_arg5 ,&ninds);  
  }                                                                                      
                                                                                         
 if (OP_diags>2) {                                                                      
    if (ninds==0)                                                                        
      printf(" kernel routine w/o indirection:  %s \n",name);                            
    else                                                                                 
      printf(" kernel routine with indirection: %s \n",name);                            
  }
  
  // initialise timers                                                            
  double cpu_t1, cpu_t2, wall_t1, wall_t2;                                        
  op_timers(&cpu_t1, &wall_t1); 
 
  if(idx0 != -1 || idx1 != -1 || idx2 != -1 || 
     idx3 != -1 || idx4 != -1 || idx5 != -1)//indirect loop
  {
      //for each indirect data set
      exchange_halo(set, arg0, acc0, idx0); 
      exchange_halo(set, arg1, acc1, idx1);
      exchange_halo(set, arg2, acc2, idx2);
      exchange_halo(set, arg3, acc3, idx3);
      exchange_halo(set, arg4, acc4, idx4);
      exchange_halo(set, arg5, acc5, idx5);
             
      //for all indirect dataset access with OP_READ
      if(acc0 == OP_READ && acc1 == OP_READ && acc2 == OP_READ && acc3 == OP_READ &&
      	 acc4 == OP_READ && acc5 == OP_READ ) exec_length = set->size;
      else  exec_length = set->size + OP_import_sets_list[set->index]->size;
  }
  else //direct loop
  {
      exec_length = set->size;
  }
                                                                                         
  // loop over set elements 
  //(1) over owned partition 
  for (int n=0; n < set->size; n++) {                                                       
    arg_set(n,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 );                     
    arg_set(n,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 );                     
    arg_set(n,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 );                     
    arg_set(n,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 );                     
    arg_set(n,arg4 ,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 );                     
    arg_set(n,arg5 ,idx5 ,map5 ,dim5 ,typ5 ,acc5 ,(T5  **)&p_arg5 );                 
    
    kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,  (T4 *)p_arg4,    
            (T5 *)p_arg5 );                                                              
  }             
  

  
  //(2) over exec halo (blank out global parameters to avoid double counting)
  for (int n=set->size; n<exec_length; n++) { 
    arg_set(n,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 );                     
    arg_set(n,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 );                     
    arg_set(n,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 );                     
    arg_set(n,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 );                     
    arg_set(n,arg4 ,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 );                     
    arg_set(n,arg5 ,idx5 ,map5 ,dim5 ,typ5 ,acc5 ,(T5  **)&p_arg5 );   
      
      // call kernel function, passing in pointers to exec halo data  
      kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,  (T4 *)p_arg4,    
            (T5 *)p_arg5 );    
  }
  
  //set dirty bit on direct/indirect datasets with access OP_INC,OP_WRITE, OP_RW
  set_dirtybit(arg0, acc0);
  set_dirtybit(arg1, acc1);
  set_dirtybit(arg2, acc2);
  set_dirtybit(arg3, acc3);
  set_dirtybit(arg4, acc4);
  set_dirtybit(arg5, acc5);

  
  //performe any global operations
  //--none
    
  //update timer record
  op_timers(&cpu_t2, &wall_t2);
  op_timing_realloc(3);
  if(OP_kernels[3].count==0)
      OP_kernels[3].name     = name;
  OP_kernels[3].count    += 1;
  OP_kernels[3].time     += wall_t2 - wall_t1;
  
}  


//                                                                                       
// op_par_loop routine for 8 arguments        - res_calc                                           
//
template < class T0, class T1, class T2, class T3, class T4,                             
           class T5, class T6, class T7,                                                 
           class D0, class D1, class D2, class D3, class D4,                             
           class D5, class D6, class D7 >                                                
void op_par_loop_res_calc(void (*kernel)( T0*, T1*, T2*, T3*, T4*,                              
                                    T5*, T6*, T7* ),                                     
  char const * name, op_set set,                                                         
  D0  arg0 ,int idx0 ,op_map map0 ,int dim0 ,const char *typ0 ,op_access acc0 ,          
  D1  arg1 ,int idx1 ,op_map map1 ,int dim1 ,const char *typ1 ,op_access acc1 ,          
  D2  arg2 ,int idx2 ,op_map map2 ,int dim2 ,const char *typ2 ,op_access acc2 ,          
  D3  arg3 ,int idx3 ,op_map map3 ,int dim3 ,const char *typ3 ,op_access acc3 ,          
  D4  arg4 ,int idx4 ,op_map map4 ,int dim4 ,const char *typ4 ,op_access acc4 ,          
  D5  arg5 ,int idx5 ,op_map map5 ,int dim5 ,const char *typ5 ,op_access acc5 ,          
  D6  arg6 ,int idx6 ,op_map map6 ,int dim6 ,const char *typ6 ,op_access acc6 ,          
  D7  arg7 ,int idx7 ,op_map map7 ,int dim7 ,const char *typ7 ,op_access acc7 ){         
                                                                                         
  char *p_arg0, *p_arg1, *p_arg2, *p_arg3, *p_arg4,                                      
       *p_arg5,*p_arg6,*p_arg7;                                                                          
                                                                                         
  int ninds=0;
  int exec_length = 0;                                                                          
                     
  
  // consistency checks                                                                  
                                                                                         
  if (OP_diags>0) {                                                                      
    arg_check(set,0 ,name,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 ,&ninds);  
    arg_check(set,1 ,name,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 ,&ninds);  
    arg_check(set,2 ,name,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 ,&ninds);  
    arg_check(set,3 ,name,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 ,&ninds);  
    arg_check(set,4 ,name,arg4 ,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 ,&ninds);  
    arg_check(set,5 ,name,arg5 ,idx5 ,map5 ,dim5 ,typ5 ,acc5 ,(T5  **)&p_arg5 ,&ninds);  
    arg_check(set,6 ,name,arg6 ,idx6 ,map6 ,dim6 ,typ6 ,acc6 ,(T6  **)&p_arg6 ,&ninds);
    arg_check(set,7 ,name,arg7 ,idx7 ,map7 ,dim7 ,typ7 ,acc7 ,(T7  **)&p_arg7 ,&ninds);  
  }                                                                                      
                                                                                         
 if (OP_diags>2) {                                                                      
    if (ninds==0)                                                                        
      printf(" kernel routine w/o indirection:  %s \n",name);                            
    else                                                                                 
      printf(" kernel routine with indirection: %s \n",name);                            
  }
  
  // initialise timers                                                            
  double cpu_t1, cpu_t2, wall_t1, wall_t2;                                        
  op_timers(&cpu_t1, &wall_t1); 
 
  if(idx0 != -1 || idx1 != -1 || idx2 != -1 || 
     idx3 != -1 || idx4 != -1 || idx5 != -1 || idx6 != -1 || idx7 != -1)//indirect loop
  {
      //for each indirect data set
      exchange_halo(set, arg0, acc0, idx0); 
      exchange_halo(set, arg1, acc1, idx1);
      exchange_halo(set, arg2, acc2, idx2);
      exchange_halo(set, arg3, acc3, idx3);
      exchange_halo(set, arg4, acc4, idx4);
      exchange_halo(set, arg5, acc5, idx5);
      exchange_halo(set, arg6, acc6, idx6);
      exchange_halo(set, arg7, acc7, idx7);
             
      //for all indirect dataset access with OP_READ
      if(acc0 == OP_READ && acc1 == OP_READ && acc2 == OP_READ && acc3 == OP_READ &&
      	 acc4 == OP_READ && acc5 == OP_READ && acc6 == OP_READ && acc7 == OP_READ )
      exec_length = set->size;
      else  exec_length = set->size + OP_import_sets_list[set->index]->size;
  }
  else //direct loop
  {
      exec_length = set->size;
  }
                                                                                         
  // loop over set elements  
  // (1) over owned partition
  for (int n=0; n < set->size; n++) {
      arg_set(n,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 );
      arg_set(n,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 );                     
      arg_set(n,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 );                     
      arg_set(n,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 );                     
      arg_set(n,arg4 ,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 );                     
      arg_set(n,arg5 ,idx5 ,map5 ,dim5 ,typ5 ,acc5 ,(T5  **)&p_arg5 );                     
      arg_set(n,arg6 ,idx6 ,map6 ,dim6 ,typ6 ,acc6 ,(T6  **)&p_arg6 );                     
      arg_set(n,arg7 ,idx7 ,map7 ,dim7 ,typ7 ,acc7 ,(T7  **)&p_arg7 );                     
                                                                                         
      // call kernel function, passing in pointers to data 
      kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,  (T4 *)p_arg4,    
            (T5 *)p_arg5,  (T6 *)p_arg6,  (T7 *)p_arg7 );  
   }
   
   //(2) over exec halo (blank out global parameters to avoid double counting)
   for (int n = set->size; n<exec_length; n++) {
      arg_set(n,arg0 ,idx0 ,map0 ,dim0 ,typ0 ,acc0 ,(T0  **)&p_arg0 );
      arg_set(n,arg1 ,idx1 ,map1 ,dim1 ,typ1 ,acc1 ,(T1  **)&p_arg1 );                     
      arg_set(n,arg2 ,idx2 ,map2 ,dim2 ,typ2 ,acc2 ,(T2  **)&p_arg2 );                     
      arg_set(n,arg3 ,idx3 ,map3 ,dim3 ,typ3 ,acc3 ,(T3  **)&p_arg3 );                     
      arg_set(n,arg4 ,idx4 ,map4 ,dim4 ,typ4 ,acc4 ,(T4  **)&p_arg4 );                     
      arg_set(n,arg5 ,idx5 ,map5 ,dim5 ,typ5 ,acc5 ,(T5  **)&p_arg5 );                     
      arg_set(n,arg6 ,idx6 ,map6 ,dim6 ,typ6 ,acc6 ,(T6  **)&p_arg6 );                     
      arg_set(n,arg7 ,idx7 ,map7 ,dim7 ,typ7 ,acc7 ,(T7  **)&p_arg7 );                    
                                                                                         
       // call kernel function, passing in pointers to exec halo data                       
       kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,  (T4 *)p_arg4,
       	   (T5 *)p_arg5,  (T6 *)p_arg6,  (T7 *)p_arg7 );  
   }
   
    //set dirty bit on direct/indirect datasets with access OP_INC,OP_WRITE, OP_RW
    set_dirtybit(arg0, acc0);
    set_dirtybit(arg1, acc1);
    set_dirtybit(arg2, acc2);
    set_dirtybit(arg3, acc3);
    set_dirtybit(arg4, acc4);
    set_dirtybit(arg5, acc5);
    set_dirtybit(arg6, acc6);
    set_dirtybit(arg7, acc7);
    
    //performe any global operations
    //--none
    
    //update timer record
    op_timers(&cpu_t2, &wall_t2);
    op_timing_realloc(2);
    if(OP_kernels[2].count==0)
    	OP_kernels[2].name     = name;
    OP_kernels[2].count    += 1;
    OP_kernels[2].time     += wall_t2 - wall_t1; 
  
}

