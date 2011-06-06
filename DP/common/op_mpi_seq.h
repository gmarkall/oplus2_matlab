//                                                                        
// header                                                                 
//                                                                        
                                                                         
                                                
                                                                          
void op_arg_set(int n, op_arg arg, char **p_arg){                         
  int n2;                                                                 
  if (arg.map==NULL)         // identity mapping, or global data          
    n2 = n;                                                               
  else                       // standard pointers                         
    n2 = arg.map->map[arg.idx+n*arg.map->dim];                            
                                                                          
  *p_arg = arg.data + n2*arg.size;                                        
}                                                                         
   

op_arg* blank_arg(op_arg *arg)
{
  
  op_arg *junck = NULL;
  if(arg->argtype == OP_ARG_GBL && arg->acc != OP_READ)//this argument is OP_GBL and OP_INC or OP_MAX/MIN
  {   
      return junck;
  }
  else 
  {
      return arg;
  }
    
}


//                                                                        
// op_par_loop routine for 2 arguments                                    
//                                                                        
                                                                          
template < class T0, class T1 >                                           
void op_par_loop(void (*kernel)( T0*, T1* ),                              
  char const * name, op_set set,                                          
  op_arg arg0, op_arg arg1 ) {                                            
                                                                          
  char *p_arg0, *p_arg1;
  int exec_length = 0;
  
  // consistency checks                                                   
                                                                          
  int ninds=0;                                                            
                                                                          
  if (OP_diags>0) {                                                       
   op_arg_check(set,0 ,arg0 ,&ninds,name);                                
   op_arg_check(set,1 ,arg1 ,&ninds,name);                                
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
  
  if(arg0.idx != -1 || arg1.idx != -1)//indirect loop
  {
      //for each indirect data set
      if(arg0.argtype == OP_ARG_DAT) exchange_halo(set, arg0);
      if(arg1.argtype == OP_ARG_DAT) exchange_halo(set, arg1);
      
      //for all indirect dataset access with OP_READ
      if(arg0.acc == OP_READ && arg1.acc == OP_READ   ) exec_length = set->size;
      else  exec_length = set->size + OP_import_sets_list[set->index]->size;
  }
  else //direct loop
  {
      exec_length = set->size;
  }     
                                                                          
  // loop over set elements                                               
  //(1) over owned partition                                                                        
  for (int n=0; n<set->size; n++) {                                       
    op_arg_set(n,arg0 ,&p_arg0 );                                         
    op_arg_set(n,arg1 ,&p_arg1 );                                         
                                                                          
    // call kernel function, passing in pointers to data 
    kernel( (T0 *)p_arg0,  (T1 *)p_arg1 );
  }
  
  //(2) over exec halo (blank out global parameters to avoid double counting)                                                                       
   for (int n=set->size; n<exec_length; n++) {                                  
    op_arg_set(n,*(blank_arg(&arg0))  ,&p_arg0 );                                         
    op_arg_set(n,*(blank_arg(&arg1))  ,&p_arg1 );                                         
                                                                          
    // call kernel function, passing in pointers to data 
    kernel( (T0 *)p_arg0,  (T1 *)p_arg1 );
  }
  
  //set dirty bit on direct/indirect datasets with access OP_INC,OP_WRITE, OP_RW
  if(arg0.argtype == OP_ARG_DAT)set_dirtybit(arg0);
  if(arg1.argtype == OP_ARG_DAT)set_dirtybit(arg1);
  
  //performe any global operations
  if(arg0.argtype == OP_ARG_GBL) 
      global_reduce(&arg0);
  if(arg1.argtype == OP_ARG_GBL) 
      global_reduce(&arg1);
  
  //update timer record
    op_timers(&cpu_t2, &wall_t2);
    
    op_timing_realloc(0);
    if(OP_kernels[0].count==0)
    	OP_kernels[0].name     = name;
    OP_kernels[0].count    += 1;
    OP_kernels[0].time     += wall_t2 - wall_t1;
}  


//                                                                        
// op_par_loop routine for 5 arguments                                    
//                                                                        
                                                                          
template < class T0, class T1, class T2, class T3,                        
           class T4 >                                                     
void op_par_loop(void (*kernel)( T0*, T1*, T2*, T3*,                      
                                 T4* ),                                   
  char const * name, op_set set,                                          
  op_arg arg0, op_arg arg1, op_arg arg2, op_arg arg3,                     
  op_arg arg4 ) {                                                         
                                                                          
  char *p_arg0, *p_arg1, *p_arg2, *p_arg3,                                
       *p_arg4;     
  int exec_length = 0;
  
  // consistency checks                                                   
                                                                          
  int ninds=0;                                                            
                                                                          
  if (OP_diags>0) {                                                       
   op_arg_check(set,0 ,arg0 ,&ninds,name);                                
   op_arg_check(set,1 ,arg1 ,&ninds,name);                                
   op_arg_check(set,2 ,arg2 ,&ninds,name);                                
   op_arg_check(set,3 ,arg3 ,&ninds,name);                                
   op_arg_check(set,4 ,arg4 ,&ninds,name);                                                            
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
  
  if(arg0.idx != -1 || arg1.idx != -1 || arg2.idx != -1 || arg3.idx != -1 ||
     arg4.idx != -1 )//indirect loop
  {
      //for each indirect data set
      if(arg0.argtype == OP_ARG_DAT) exchange_halo(set, arg0);
      if(arg1.argtype == OP_ARG_DAT) exchange_halo(set, arg1);
      if(arg2.argtype == OP_ARG_DAT) exchange_halo(set, arg2);
      if(arg3.argtype == OP_ARG_DAT) exchange_halo(set, arg3);
      if(arg4.argtype == OP_ARG_DAT) exchange_halo(set, arg4);

      //for all indirect dataset access with OP_READ
      if(arg0.acc == OP_READ && arg1.acc == OP_READ && arg2.acc == OP_READ && arg3.acc == OP_READ &&
      	  arg4.acc == OP_READ )
      exec_length = set->size;
      else  exec_length = set->size + OP_import_sets_list[set->index]->size;
  }
  else //direct loop
  {
      exec_length = set->size;
  }
                                                                               
  // loop over set elements                                               
  //(1) over owned partition                                                                        
  for (int n=0; n<set->size; n++) {                                       
    op_arg_set(n,arg0 ,&p_arg0 );                                         
    op_arg_set(n,arg1 ,&p_arg1 );                                         
    op_arg_set(n,arg2 ,&p_arg2 );                                         
    op_arg_set(n,arg3 ,&p_arg3 );                                         
    op_arg_set(n,arg4 ,&p_arg4 );                                                                               
                                                                          
    // call kernel function, passing in pointers to data 
    kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,    
            (T4 *)p_arg4 );  
  }
  
  

  //(2) over exec halo (blank out global parameters to avoid double counting)                                                                       
   for (int n=set->size; n<exec_length; n++) {                                  
       op_arg_set(n,*(blank_arg(&arg0)) ,&p_arg0 );
       op_arg_set(n,*(blank_arg(&arg1)) ,&p_arg1 );
       op_arg_set(n,*(blank_arg(&arg2)) ,&p_arg2 );
       op_arg_set(n,*(blank_arg(&arg3)) ,&p_arg3 );
       op_arg_set(n,*(blank_arg(&arg4)) ,&p_arg4 );
                                                                          
    // call kernel function, passing in pointers to data 
    kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,    
            (T4 *)p_arg4 );  
  }
  
  //set dirty bit on direct/indirect datasets with access OP_INC,OP_WRITE, OP_RW
  if(arg0.argtype == OP_ARG_DAT)set_dirtybit(arg0);
  if(arg1.argtype == OP_ARG_DAT)set_dirtybit(arg1);
  if(arg2.argtype == OP_ARG_DAT)set_dirtybit(arg2);
  if(arg3.argtype == OP_ARG_DAT)set_dirtybit(arg3);
  if(arg4.argtype == OP_ARG_DAT)set_dirtybit(arg4);
  
  //performe any global operations
  if(arg0.argtype == OP_ARG_GBL) 
      global_reduce(&arg0);
  if(arg1.argtype == OP_ARG_GBL) 
      global_reduce(&arg1);;
  if(arg2.argtype == OP_ARG_GBL) 
      global_reduce(&arg2);
  if(arg3.argtype == OP_ARG_GBL) 
      global_reduce(&arg3);
  if(arg4.argtype == OP_ARG_GBL) 
      global_reduce(&arg4);
  
  //update timer record
  op_timers(&cpu_t2, &wall_t2);
  
  op_timing_realloc(4);
  if(OP_kernels[4].count==0)
      OP_kernels[4].name     = name;
  OP_kernels[4].count    += 1;
  OP_kernels[4].time     += wall_t2 - wall_t1;
}  



//                                                                        
// op_par_loop routine for 6 arguments                                    
//                                                                        
                                                                          
template < class T0, class T1, class T2, class T3,                        
           class T4, class T5 >                                           
void op_par_loop(void (*kernel)( T0*, T1*, T2*, T3*,                      
                                 T4*, T5* ),                              
  char const * name, op_set set,                                          
  op_arg arg0, op_arg arg1, op_arg arg2, op_arg arg3,                     
  op_arg arg4, op_arg arg5 ) {                                            
                                                                          
  char *p_arg0, *p_arg1, *p_arg2, *p_arg3,                                
       *p_arg4, *p_arg5;    
  int exec_length = 0;
  
  // consistency checks                                                   
                                                                          
  int ninds=0;                                                            
                                                                          
  if (OP_diags>0) {                                                       
   op_arg_check(set,0 ,arg0 ,&ninds,name);                                
   op_arg_check(set,1 ,arg1 ,&ninds,name);                                
   op_arg_check(set,2 ,arg2 ,&ninds,name);                                
   op_arg_check(set,3 ,arg3 ,&ninds,name);                                
   op_arg_check(set,4 ,arg4 ,&ninds,name);                                
   op_arg_check(set,5 ,arg5 ,&ninds,name);                                
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
  
  if(arg0.idx != -1 || arg1.idx != -1 || arg2.idx != -1 || arg3.idx != -1 ||
     arg4.idx != -1 || arg5.idx != -1)//indirect loop
  {
      //for each indirect data set
      if(arg0.argtype == OP_ARG_DAT)exchange_halo(set, arg0);
      if(arg1.argtype == OP_ARG_DAT)exchange_halo(set, arg1);
      if(arg2.argtype == OP_ARG_DAT)exchange_halo(set, arg2);
      if(arg3.argtype == OP_ARG_DAT)exchange_halo(set, arg3);
      if(arg4.argtype == OP_ARG_DAT)exchange_halo(set, arg4);
      if(arg5.argtype == OP_ARG_DAT)exchange_halo(set, arg5);

      //for all indirect dataset access with OP_READ
      if(arg0.acc == OP_READ && arg1.acc == OP_READ && arg2.acc == OP_READ && arg3.acc == OP_READ &&
      	  arg4.acc == OP_READ && arg5.acc == OP_READ)
      exec_length = set->size;
      else  exec_length = set->size + OP_import_sets_list[set->index]->size;
  }
  else //direct loop
  {
      exec_length = set->size;
  }
                                                                               
  // loop over set elements                                               
  //(1) over owned partition                                                                        
  for (int n=0; n<set->size; n++) {                                       
    op_arg_set(n,arg0 ,&p_arg0 );                                         
    op_arg_set(n,arg1 ,&p_arg1 );                                         
    op_arg_set(n,arg2 ,&p_arg2 );                                         
    op_arg_set(n,arg3 ,&p_arg3 );                                         
    op_arg_set(n,arg4 ,&p_arg4 );                                         
    op_arg_set(n,arg5 ,&p_arg5 );                                         
                                                                          
    // call kernel function, passing in pointers to data 
    kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,    
            (T4 *)p_arg4,  (T5 *)p_arg5 );  
  }
  
  //(2) over exec halo (blank out global parameters to avoid double counting)                                                                       
   for (int n=set->size; n<exec_length; n++) {                                  
    op_arg_set(n,*(blank_arg(&arg0)) ,&p_arg0 );                                         
    op_arg_set(n,*(blank_arg(&arg1)) ,&p_arg1 );                                         
    op_arg_set(n,*(blank_arg(&arg2)) ,&p_arg2 );                                         
    op_arg_set(n,*(blank_arg(&arg3)) ,&p_arg3 );                                         
    op_arg_set(n,*(blank_arg(&arg4)) ,&p_arg4 );                                         
    op_arg_set(n,*(blank_arg(&arg5)) ,&p_arg5 );                                         
                                                                          
    // call kernel function, passing in pointers to data 
    kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,    
            (T4 *)p_arg4,  (T5 *)p_arg5 );  
  }
  
  //set dirty bit on direct/indirect datasets with access OP_INC,OP_WRITE, OP_RW
  if(arg0.argtype == OP_ARG_DAT)set_dirtybit(arg0);
  if(arg1.argtype == OP_ARG_DAT)set_dirtybit(arg1);
  if(arg2.argtype == OP_ARG_DAT)set_dirtybit(arg2);
  if(arg3.argtype == OP_ARG_DAT)set_dirtybit(arg3);
  if(arg4.argtype == OP_ARG_DAT)set_dirtybit(arg4);
  if(arg5.argtype == OP_ARG_DAT)set_dirtybit(arg5);
  
  //performe any global operations
  if(arg0.argtype == OP_ARG_GBL) 
      global_reduce(&arg0);
  if(arg1.argtype == OP_ARG_GBL) 
      global_reduce(&arg1);;
  if(arg2.argtype == OP_ARG_GBL) 
      global_reduce(&arg2);
  if(arg3.argtype == OP_ARG_GBL) 
      global_reduce(&arg3);
  if(arg4.argtype == OP_ARG_GBL) 
      global_reduce(&arg4);
  if(arg5.argtype == OP_ARG_GBL) 
      global_reduce(&arg5);
  
  //update timer record
    op_timers(&cpu_t2, &wall_t2);
    
    op_timing_realloc(1);
    if(OP_kernels[1].count==0)
    	OP_kernels[1].name     = name;
    OP_kernels[1].count    += 1;
    OP_kernels[1].time     += wall_t2 - wall_t1;
}  




//                                                                        
// op_par_loop routine for 8 arguments                                    
//                                                                        
                                                                          
template < class T0, class T1, class T2, class T3,                        
           class T4, class T5, class T6, class T7 >                       
void op_par_loop(void (*kernel)( T0*, T1*, T2*, T3*,                      
                                 T4*, T5*, T6*, T7* ),                    
  char const * name, op_set set,                                          
  op_arg arg0, op_arg arg1, op_arg arg2, op_arg arg3,                     
  op_arg arg4, op_arg arg5, op_arg arg6, op_arg arg7 ) {                  
                                                                          
  char *p_arg0, *p_arg1, *p_arg2, *p_arg3,                                
       *p_arg4, *p_arg5, *p_arg6, *p_arg7;          
  int exec_length = 0;
  
  // consistency checks                                                   
                                                                          
  int ninds=0;                                                            
                                                                          
  if (OP_diags>0) {                                                       
   op_arg_check(set,0 ,arg0 ,&ninds,name);                                
   op_arg_check(set,1 ,arg1 ,&ninds,name);                                
   op_arg_check(set,2 ,arg2 ,&ninds,name);                                
   op_arg_check(set,3 ,arg3 ,&ninds,name);                                
   op_arg_check(set,4 ,arg4 ,&ninds,name);                                
   op_arg_check(set,5 ,arg5 ,&ninds,name);                                
   op_arg_check(set,6 ,arg6 ,&ninds,name);                                
   op_arg_check(set,7 ,arg7 ,&ninds,name);                              
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
  
  if(arg0.idx != -1 || arg1.idx != -1 || arg2.idx != -1 || arg3.idx != -1 ||
     arg4.idx != -1 || arg5.idx != -1 || arg6.idx != -1 || arg7.idx != -1)//indirect loop
  {
      //for each indirect data set
      if(arg0.argtype == OP_ARG_DAT)exchange_halo(set, arg0);
      if(arg1.argtype == OP_ARG_DAT)exchange_halo(set, arg1);
      if(arg2.argtype == OP_ARG_DAT)exchange_halo(set, arg2);
      if(arg3.argtype == OP_ARG_DAT)exchange_halo(set, arg3);
      if(arg4.argtype == OP_ARG_DAT)exchange_halo(set, arg4);
      if(arg5.argtype == OP_ARG_DAT)exchange_halo(set, arg5);
      if(arg6.argtype == OP_ARG_DAT)exchange_halo(set, arg6);
      if(arg7.argtype == OP_ARG_DAT)exchange_halo(set, arg7);
      
      //for all indirect dataset access with OP_READ
      if(arg0.acc == OP_READ && arg1.acc == OP_READ && arg2.acc == OP_READ && arg3.acc == OP_READ &&
      	  arg4.acc == OP_READ && arg5.acc == OP_READ && arg6.acc == OP_READ && arg7.acc == OP_READ)
      exec_length = set->size;
      else  exec_length = set->size + OP_import_sets_list[set->index]->size;
  }
  else //direct loop
  {
      exec_length = set->size;
  }
                                                                               
  // loop over set elements                                               
  //(1) over owned partition                                                                        
  for (int n=0; n<set->size; n++) {                                       
    op_arg_set(n,arg0 ,&p_arg0 );                                         
    op_arg_set(n,arg1 ,&p_arg1 );                                         
    op_arg_set(n,arg2 ,&p_arg2 );                                         
    op_arg_set(n,arg3 ,&p_arg3 );                                         
    op_arg_set(n,arg4 ,&p_arg4 );                                         
    op_arg_set(n,arg5 ,&p_arg5 );                                         
    op_arg_set(n,arg6 ,&p_arg6 );                                         
    op_arg_set(n,arg7 ,&p_arg7 );                                             
                                                                          
    // call kernel function, passing in pointers to data 
    kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,    
            (T4 *)p_arg4,  (T5 *)p_arg5,  (T6 *)p_arg6,  (T7 *)p_arg7 );  
  }
  
  //(2) over exec halo (blank out global parameters to avoid double counting)                                                                       
   for (int n=set->size; n<exec_length; n++) {                                  
    op_arg_set(n,*(blank_arg(&arg0)) ,&p_arg0 );                                         
    op_arg_set(n,*(blank_arg(&arg1)) ,&p_arg1 );                                         
    op_arg_set(n,*(blank_arg(&arg2)) ,&p_arg2 );                                         
    op_arg_set(n,*(blank_arg(&arg3)) ,&p_arg3 );                                         
    op_arg_set(n,*(blank_arg(&arg4)) ,&p_arg4 );                                         
    op_arg_set(n,*(blank_arg(&arg5)) ,&p_arg5 );                                            
    op_arg_set(n,*(blank_arg(&arg6)) ,&p_arg6 );                                         
    op_arg_set(n,*(blank_arg(&arg7)) ,&p_arg7 );                                             
                                                                          
    // call kernel function, passing in pointers to data 
    kernel( (T0 *)p_arg0,  (T1 *)p_arg1,  (T2 *)p_arg2,  (T3 *)p_arg3,    
            (T4 *)p_arg4,  (T5 *)p_arg5,  (T6 *)p_arg6,  (T7 *)p_arg7 );  
  }
  
  //set dirty bit on direct/indirect datasets with access OP_INC,OP_WRITE, OP_RW
  if(arg0.argtype == OP_ARG_DAT)set_dirtybit(arg0);
  if(arg1.argtype == OP_ARG_DAT)set_dirtybit(arg1);
  if(arg2.argtype == OP_ARG_DAT)set_dirtybit(arg2);
  if(arg3.argtype == OP_ARG_DAT)set_dirtybit(arg3);
  if(arg4.argtype == OP_ARG_DAT)set_dirtybit(arg4);
  if(arg5.argtype == OP_ARG_DAT)set_dirtybit(arg5);
  if(arg6.argtype == OP_ARG_DAT)set_dirtybit(arg6);
  if(arg7.argtype == OP_ARG_DAT)set_dirtybit(arg7);
  
  //performe any global operations
  if(arg0.argtype == OP_ARG_GBL) 
      global_reduce(&arg0);
  if(arg1.argtype == OP_ARG_GBL) 
      global_reduce(&arg1);;
  if(arg2.argtype == OP_ARG_GBL) 
      global_reduce(&arg2);
  if(arg3.argtype == OP_ARG_GBL) 
      global_reduce(&arg3);
  if(arg4.argtype == OP_ARG_GBL) 
      global_reduce(&arg4);
  if(arg5.argtype == OP_ARG_GBL) 
      global_reduce(&arg5);
  if(arg6.argtype == OP_ARG_GBL) 
      global_reduce(&arg6);
  if(arg7.argtype == OP_ARG_GBL) 
      global_reduce(&arg7);
  
  //update timer record
    op_timers(&cpu_t2, &wall_t2);    
    
    
    op_timing_realloc(2);
    if(OP_kernels[2].count==0)
    	OP_kernels[2].name     = name;
    OP_kernels[2].count    += 1;
    OP_kernels[2].time     += wall_t2 - wall_t1;
}  
