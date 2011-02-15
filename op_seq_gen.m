

function op_seq_gen()

%
% this sets the max number of arguments in op_par_loop
%

maxargs = 20;

%
% first the top bit with headers and function prototypes
%

file = strvcat(...',...
'//',...
'// headers',...
'//',...
' ',...
'#include <stdlib.h>',...
'#include <stdio.h>',...
'#include <string.h>',...
'#include <math.h>',...
'#include "op_datatypes.h"',...
' ',...
'//',...
'// op routine declarations',...
'//',...
' ',...
'void op_init(int, char **);',...
' ',...
'void op_decl_set(int, op_set &, char const *);',...
' ',...
'void op_decl_ptr(op_set, op_set, int, int *, op_ptr &, char const *);',...
' ',...
'void op_decl_dat(op_set, int, op_datatype, double *, op_dat &, char const *);',...
'void op_decl_dat(op_set, int, op_datatype, float  *, op_dat &, char const *);',...
'void op_decl_dat(op_set, int, op_datatype, int    *, op_dat &, char const *);',...
' ',...
'void op_decl_const(int, op_datatype, double *, char const *);',...
'void op_decl_const(int, op_datatype, float  *, char const *);',...
'void op_decl_const(int, op_datatype, int    *, char const *);',...
' ',...
'void op_fetch_data(op_dat);',...
' ',...
'void op_diagnostic_output();',...
' ',' ',...
'template < class T, class D >',...
'void arg_check(op_set set, int m, const char *name, D * arg, int idx, op_ptr ptr,',...
'               int dim, op_datatype typ, op_access acc, T * p_arg, int *ninds){',...
'  if (idx != -1) {',...
'    printf("error: arg %d in kernel \"%s\"\n",m,name);',...
'    printf("invalid index, should be -1 for constant arg\n");',...
'    exit(1);',...
'  }',...
'  if (type_error(arg,typ)) {',...
'    printf("error: arg %d in kernel \"%s\"\n",m,name);',...
'    printf("OP_datatype does not match type of input dataset\n");',...
'    exit(1);',...
'  }',...
'  if (type_error(p_arg,typ)) {',...
'    printf("error: arg %d in kernel \"%s\"\n",m,name);',...
'    printf("OP_datatype does not match type of function argument\n");',...
'    exit(1);',...
'  }',...
'  if (dim <= 0) {',...
'    printf("error: arg %d in kernel \"%s\"\n",m,name);',...
'    printf("dimension must be strictly positive\n");',...
'    exit(1);',...
'  }',...
'}',...
' ',' ',...
'template < class T >',...
'void arg_check(op_set set, int m, const char *name, op_dat arg, int idx, op_ptr ptr,',...
'               int dim, op_datatype typ, op_access acc, T * p_arg, int *ninds){',...
'  if (idx>=0) (*ninds)++;',...
' ',...
'  if (ptr.ptr == NULL) {',...
'    if (idx != -1) {',...
'      printf("error: arg %d in kernel \"%s\"\n",m,name);',...
'      printf("invalid index, should be -1 for identity mapping\n");',...
'      exit(1);',...
'    }',...
'  }',...                                                                                  
'  else {',...
'    if (set.index != ptr.from.index || arg.set.index != ptr.to.index) {',...
'      printf("error: arg %d in kernel \"%s\"\n",m,name);',...
'      printf("invalid pointer\n");',...
'      exit(1);',...
'    }',...
'    if (idx < 0 || idx >= ptr.dim) {',...
'      printf("error: arg %d in kernel \"%s\"\n",m,name);',...
'      printf("invalid pointer index\n");',...
'      exit(1);',...
'    }',...
'  }',...
'  if (arg.type != typ) {',...
'    printf("error: arg %d in kernel \"%s\"\n",m,name);',...
'    printf("OP_datatype does not match type of input dataset\n");',...
'    exit(1);',...
'  }',...
'  if (type_error(p_arg,typ)) {',...
'    printf("error: arg %d in kernel \"%s\"\n",m,name);',...
'    printf("OP_datatype does not match type of function argument\n");',...
'    exit(1);',...
'  }',...
'  if (arg.dim != dim) {',...
'    printf("error: arg %d in kernel \"%s\"\n",m,name);',...
'    printf("dimension does not match input dataset\n");',...
'    exit(1);',...
'  }',...
'}',...
' ',' ',...
'template < class T >',...
'void arg_set(int n,T *arg,int idx,op_ptr ptr,int dim,',...
'             op_datatype typ,op_access acc,char **p_arg){',...
'  *p_arg = (char *) arg;',...
'}',...
' ',' ',...
'void arg_set(int n,op_dat arg,int idx,op_ptr ptr,int dim,',...
'             op_datatype typ,op_access acc,char **p_arg){',...
'  int n2;',...
'  if (ptr.dim==0)                   // identity mapping',...
'    n2 = n;',...
'  else                              // standard pointers',...
'    n2 = ptr.ptr[idx+n*ptr.dim];',...
' ',...
'  *p_arg = arg.dat + n2*arg.size;',...
'}',...
' ');

%
% now for op_par_loop defns
%

for nargs = 1:maxargs

  file = strvcat(file,' ', ...
    '// ',...
   ['// op_par_loop routine for ' num2str(nargs) ' arguments '],...
    '// ',' ');

  n_per_line = 5;

  line = 'template < ';
  for n = 1:nargs
    line = [ line 'class T' num2str(n-1) ','];
    if (mod(n,n_per_line)==0 || n==nargs)
      file = strvcat(file,line);
      line = '           ';
    elseif (n<=10)
      line = [line ' '];
    end
  end
  for n = 1:nargs
    line = [ line 'class D' num2str(n-1) ','];
    if (n==nargs)
      line = [line(1:end-1) ' >'];
    end
    if (mod(n,n_per_line)==0 || n==nargs)
      file = strvcat(file,line);
      line = '           ';
    elseif (n<=10)
      line = [line ' '];
    end
  end

  line = ['void op_par_loop_' num2str(nargs) '(void (*kernel)( '];
  for n = 1:nargs
    line = [ line 'T' num2str(n-1) '*,'];
    if (n==nargs) 
      line = [line(1:end-1) ' ),'];
    end
    if (mod(n,n_per_line)==0 || n==nargs)
      file = strvcat(file,line);
      line = '                                    ';
    elseif (n<=10)
      line = [line ' '];
    end
  end

  file = strvcat(file,'  char const * name, op_set set,');

  for n = 1:nargs
    if (n<=10)
      line = '  D0  arg0, int idx0, op_ptr ptr0, int dim0, op_datatype typ0, op_access acc0,';
    else
      line = '  D0 arg0,int idx0,op_ptr ptr0,int dim0,op_datatype typ0,op_access acc0,';
    end
    line = regexprep(line,'0',num2str(n-1));
    if (n==nargs)
      line = [line(1:end-1) '){'];
    end
    file = strvcat(file,line);
  end

  file = strvcat(file,' ');

  for n = 1:nargs
    line = [ '  char *p_arg' num2str(n-1) ';' ];
    file = strvcat(file,line);
  end

%
% diagnostics
%

  file = strvcat(file,' ',...
                      '  int ninds=0;',' ',...
                      '  // consistency checks',' ',...
                      '  if (OP_DIAGS>0) {');

  for n = 1:nargs
    line = '    arg_check(set, 0,name, arg0,idx0,ptr0,dim0,typ0,acc0,(T0 *)p_arg0,&ninds);';
    line = regexprep(line,'0',num2str(n-1));
    file = strvcat(file,line);
  end

  file = strvcat(file,'  }',' ',...
                      '  if (OP_DIAGS>1) {',...
                      '    if (ninds==0)',...
                      '      printf(" kernel routine w/o indirection:  %s \n",name);',...
                      '    else',...
                      '      printf(" kernel routine with indirection: %s \n",name);',...
                      '  }',' ');

%
% main loop
%

  file = strvcat(file,' ',...
                      '  // loop over set elements',' ',...
                      '  for (int n=0; n<set.size; n++) {');
  for n = 1:nargs
    line = '    arg_set(n,arg0,idx0,ptr0,dim0,typ0,acc0, &p_arg0);';
    line = regexprep(line,'0',num2str(n-1));
    file = strvcat(file,line);
  end

%
% call to user's kernel
%

  n_per_line = 5;

  file = strvcat(file,' ','    // call kernel function, passing in pointers to data',' ');

  line = ['    kernel( '];
  for n = 1:nargs
    line = [ line '(T' num2str(n-1) ' *)p_arg'  num2str(n-1) ','];
    if (n==nargs) 
      line = [line(1:end-1) ' );'];
    end
    if (mod(n,n_per_line)==0 || n==nargs)
      file = strvcat(file,line);
      line = '            ';
    elseif (n<=10)
      line = [line '  '];
    end
  end

  file = strvcat(file,'  }','}');

end


%
% print out into file
%


fid = fopen('op_seq.h','wt');
fprintf(fid,'// \n// auto-generated by op_seq_gen.m on %s \n//\n\n',datestr(now));
for n=1:size(file,1)
  fprintf(fid,'%s\n',file(n,:));
end
fclose(fid);

