%
% Source code transformation tool
%
% This tool parses the user's original source code
% to produce the CUDA code which will execute the
% user's kernel functions.
%
% The main deficiency of this current implementation
% is that it does not handle cases in which multiple 
% arguments involve the same underlying dataset.
%
% This prototype is written in MATLAB but a future
% version may use Python to avoid licensing costs.
% Alternatively, the MATLAB processor can be 
% "compiled" to produce a standalone version which
% can be freely distributed.
%
%
% usage: op2('filename')
%
% This takes as input 
%
% filename.cpp
%
% and produces as output
%
% filename_op.cpp
% filename_kernels.cu
%
% and one or more files of the form
%
% xxx_kernel.cu
%
% where xxx corresponds to the name of one of the
% kernel functions in filename.cpp
%

function op2(filename)

global dims idxs typs indtyps inddims

global OP_FLOAT OP_DOUBLE OP_INT
global OP_READ OP_WRITE OP_RW OP_INC

%
% declare constants
%

OP_FLOAT  = 1;
OP_DOUBLE = 2;
OP_INT    = 3;

OP_READ  = 1;
OP_WRITE = 2;
OP_RW    = 3;
OP_INC   = 4;

%
% read in source file and strip out white space
%

src_file = fileread([filename '.cpp']);
src_file = regexprep(src_file,'\s','');

nkernels = 0;

while (~isempty(strfind(src_file,'op_par_loop_')))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% parse file for next op_par_loop
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  loc  = min(strfind(src_file,'op_par_loop_'));
  src_file = src_file(loc+12:end);

  [num,  src_file] = strtok(src_file,'(');
  nargs = str2num(num);
  [src_args, src_file] = strtok(src_file,')');
  src_args = src_args(2:end);

  loc = [0 strfind(src_args,',') length(src_args)+1];

  na = length(loc)-1;

  if( na ~= 3+6*nargs)
    disp(sprintf('wrong number of arguments: expected %d, found %d',...
         3+6*nargs, na));
    error('aborting')
  end

  for n = 1:na
    C{n} = src_args(loc(n)+1:loc(n+1)-1);
  end

  nkernels = nkernels + 1;
  kernel{nkernels} = C{1};
  kernel_args(nkernels) = nargs;

%
% process parameters
%

  fn_name = C{1};

  for m = 1:3
    idxs(m) = str2num(C{-1+6*m});
    dims(m) = str2num(C{1+6*m});

    switch C{2+6*m}
      case 'OP_FLOAT'
        typs(m) = OP_FLOAT;
      case 'OP_DOUBLE'
        typs(m) = OP_DOUBLE;
      case 'OP_INT'
        typs(m) = OP_INT;
      otherwise
        disp('unknown data type')
    end

    switch C{3+6*m}
      case 'OP_READ'
        accs(m) = OP_READ;
      case 'OP_WRITE'
        accs(m) = OP_WRITE;
      case 'OP_RW'
        accs(m) = OP_RW;
      case 'OP_INC'
        accs(m) = OP_INC;
      otherwise
        disp('unknown access type')
    end
  end

%  disp(sprintf('name =  %s',fn_name))
%  disp(strcat('idxs =  ',sprintf('  %d',idxs)))
%  disp(strcat('dims =  ',sprintf('  %d',dims)))
%  disp(strcat('typs =  ',sprintf('  %d',typs)))
%  disp(strcat('accs =  ',sprintf('  %d',accs)))
%  disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% identify indirect datasets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ninds = 0;
  invinds  = zeros(1,nargs);
  inds     = zeros(1,nargs);
  indtyps  = ones(1,nargs);
  inddims  = zeros(1,nargs);
  indaccs  = zeros(1,nargs);

  for type = [OP_DOUBLE OP_FLOAT OP_INT]
    j = find( idxs>=0 & typs==type );
    while (~isempty(j))
      match = strcmp(C(-2+6*j(1)),C(-2+6*j)) & (accs(j(1))==accs(j));
      ninds = ninds + 1; 
      indtyps(ninds) = typs(j(1));
      inddims(ninds) = dims(j(1));
      indaccs(ninds) = accs(j(1));
      inds(j(find(match))) = ninds;
      invinds(ninds) = j(1);
      j = j(find(~match));
    end
  end

  % disp([sprintf(' indirect datasets = %d, mapping: ',ninds)  num2str(inds)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% create new file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  file = strvcat('// header files          ',...
                 '                         ',...
                 '#include <stdlib.h>      ',...
                 '#include <stdio.h>       ',...
                 '#include <string.h>      ',...
                 '#include <math.h>        ',...
                 '#include <cutil_inline.h>',...
                 '#include "op_datatypes.h"',...
                 '                         ',...
                 '                         ',...
                 'extern "C"               ',...
                 'op_plan * plan(char const *, op_set, int, op_dat *, int *,',...
                 '  op_ptr *, int *, op_datatype *, op_access *, int, int *);',...
                 '                         ',...
                 '                         ',...
                 '// user function         ',...
                 '                         ',...
                 '__device__               ');
  file = strvcat(file,['#include "' fn_name '.h"'],' ',' ');
  file = strvcat(file,'// CUDA kernel function',' ');

  if (max(idxs)>=0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% standard version 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  file = strvcat(file,strcat('__global__ void op_cuda_',fn_name,'('));

  for m = 1:ninds
    line = '  INDTYP *ind_ARG, int *ind_ARG_ptr, int *ind_ARG_sizes, int *ind_ARG_offset,';
    file = strvcat(file,rep(line,m));
  end

  for m = 1:nargs
    if (idxs(m)<0)
      line = '  TYP *ARG,';
    else
      line = '  int   *ARG_ptr,';
    end
    file = strvcat(file,rep(line,m));
  end

  file = strvcat(file,'  int    block_offset,',...
                      '  int   *blkmap,      ',...
                      '  int   *offset,      ',...
                      '  int   *nelems,      ',...
                      '  int   *ncolors,     ',...
                      '  int   *colors) {    ',' ');

  for m = 1:ninds
    line = '  INDTYP *ind_ARG_s;';
    file = strvcat(file,rep(line,m));
  end

  for m = 1:nargs
    line = '  TYP  ARG_l[DIM];';
    file = strvcat(file,rep(line,m));
  end

  file = strvcat(file,' ',...
                      '  extern __shared__ int shared[];',' ',...
                      '  // get block info',' ',...
                      '  int blockId = blkmap[blockIdx.x + block_offset];',' ');

  file = strvcat(file,'  // get sizes and shift pointers and direct-mapped data',' ');
  for m = 1:ninds
    line = '  int ind_ARG_size = ind_ARG_sizes[blockId];';
    file = strvcat(file,rep(line,m));
  end
  file = strvcat(file,' ');
  for m = 1:ninds
    line = '  ind_ARG_ptr += ind_ARG_offset[blockId];';
    file = strvcat(file,rep(line,m));
  end
  file = strvcat(file,' ');
  for m = 1:nargs
    if (idxs(m)<0)
      line = '  ARG         += offset[blockId]*DIM;';
    else
      line = '  ARG_ptr     += offset[blockId];';
    end
    file = strvcat(file,rep(line,m));
  end

  file = strvcat(file,'  colors       += offset[blockId];',' ');

  file = strvcat(file,'  // set shared memory pointers',' ');
  for m = 1:ninds
    if (m==1)
      line = '  ind_ARG_s = (INDTYP *) shared;';
      file = strvcat(file,rep(line,m));
    else
      line1 = '  ind_ARG_s = (INDTYP *)';
      line2 = ' &ind_ARG_s[ind_ARG_size*INDDIM];';
      file = strvcat(file,[rep(line1,m) rep(line2,m-1)]);
    end
  end

  file = strvcat(file,' ','  // copy indirect datasets into shared memory or zero increment',' ');
  for m = 1:ninds
    line = '  for (int n=threadIdx.x; n<INDARG_size; n+=blockDim.x)';
    file = strvcat(file,rep(line,m));
    line = '    for (int d=0; d<INDDIM; d++)';
    file = strvcat(file,rep(line,m));
    if(indaccs(m)==OP_READ | indaccs(m)==OP_RW)
      line = '      INDARG_s[d+n*INDDIM] = INDARG[d+INDARG_ptr[n]*INDDIM];';
    elseif(indaccs(m)==OP_INC)
      line = '      INDARG_s[d+n*INDDIM] = 0;';
    end
    file = strvcat(file,rep(line,m),' ');
  end
  file = strvcat(file,'  __syncthreads();',' ');

  file = strvcat(file,'  // process set elements                     ',' ', ...
                      '  int nelems2 = blockDim.x*(1+(nelems[blockId]-1)/blockDim.x);',' ',...
                      '  for (int n=threadIdx.x; n<nelems2; n+=blockDim.x) {         ',    ...
                      '    int col2 = ncolors[blockId];                              ',' ',...
                      '    if (n<nelems[blockId]) {                                  ',' ',...
                      '      // initialise local variables                           ',' ');

  for m = 1:nargs
    if (accs(m)==OP_READ | accs(m)==OP_RW)
      line = '      for (int d=0; d<DIM; d++)';
      file = strvcat(file,rep(line,m));
      if (idxs(m)<0)
        line = '        ARG_l[d] = ARG[d+n*DIM];';
      else
        line = sprintf('        ARG_l[d] = ind_arg%d_s[d+ARG_ptr[n]*DIM];',inds(m)-1);
      end
      file = strvcat(file,rep(line,m));

    elseif (accs(m)==OP_INC)
      line = '      for (int d=0; d<DIM; d++)';
      file = strvcat(file,rep(line,m));
      line = '        ARG_l[d] = 0;';
      file = strvcat(file,rep(line,m));
    end
  end

  file = strvcat(file,' ','      // user-supplied kernel call',' ');

  line = regexprep('      FN(','FN',fn_name);
  for m = 1:nargs
    line = strcat(line,rep('ARG_l,',m));
  end
  line = strcat(line(1:end-1),');');
  file = strvcat(file,line,' ');

  file = strvcat(file,'       col2 = colors[n];     ',' ',...
                      '    }                        ',' ',...
                      '    // store local variables ',' ');

  for m = 1:nargs
    if (accs(m)==OP_WRITE | accs(m)==OP_RW)
      line = '    for (int d=0; d<DIM; d++)';
      file = strvcat(file,rep(line,m));
      if (idxs(m)<0)
        line = '      ARG[d+n*DIM] = ARG_l[d];';
      else
        line = sprintf('      ind_arg%d_s[d+ARG_ptr[n]*DIM] = ARG_l[d];',inds(m)-1);
      end
      file = strvcat(file,rep(line,m));

    elseif (accs(m)==OP_INC)
      if (idxs(m)<0)
        line = '    for (int d=0; d<DIM; d++)';
        file = strvcat(file,rep(line,m));
        line = '      ARG[d+n*DIM] += ARG_l[d];';
        file = strvcat(file,rep(line,m));
      else
        file = strvcat(file, ...
         '    for (int col=0; col<ncolors[blockId]; col++) {', ...
         '      if (col2==col) ');
        line = '        for (int d=0; d<DIM; d++)';
        file = strvcat(file,rep(line,m));
        line = sprintf('          ind_arg%d_s[d+ARG_ptr[n]*DIM] += ARG_l[d];',inds(m)-1);
        file = strvcat(file,rep(line,m), ...
         '      __syncthreads();','    }',' ');
      end
    end
  end

  file = strvcat(file,'  }',' ','  // apply pointered write/increment',' ');
  for m = 1:ninds
    if(indaccs(m)==OP_WRITE | indaccs(m)==OP_RW)
      line = '  for (int n=threadIdx.x; n<INDARG_size; n+=blockDim.x)';
      file = strvcat(file,rep(line,m));
      line = '    for (int d=0; d<INDDIM; d++)';
      file = strvcat(file,rep(line,m));
      line = '      INDARG[d+INDARG_ptr[n]*INDDIM] = INDARG_s[d+n*DIM];';
      file = strvcat(file,rep(line,m),' ');
    elseif(indaccs(m)==OP_INC)
      line = '  for (int n=threadIdx.x; n<INDARG_size; n+=blockDim.x)';
      file = strvcat(file,rep(line,m));
      line = '    for (int d=0; d<INDDIM; d++)';
      file = strvcat(file,rep(line,m));
      line = '      INDARG[d+INDARG_ptr[n]*INDDIM] += INDARG_s[d+n*DIM];';
      file = strvcat(file,rep(line,m),' ');
    end
  end

  file = strvcat(file,'}')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% add stub function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  file = strvcat(file,' ',' ','// host stub function            ',' ','extern "C"');
  file = strvcat(file,['void op_par_loop_' fn_name '(char const * name, op_set set,']);
  for m = 1:nargs
    line = ['  op_dat arg%d, int idx%d, op_ptr ptr%d, int dim%d,' ...
            ' op_datatype typ%d, op_access acc%d'];
    line = sprintf(line, m-1,m-1,m-1,m-1,m-1,m-1);
    if (m<nargs)
      file = strvcat(file,[line ',']);
    else
      file = strvcat(file,[line '){']);
    end
  end

  file = strvcat(file,' ',sprintf('  int         nargs = %d, ninds = %d;',nargs,ninds),' ');

  for l=1:6
    if (l==1)
      word = 'arg'
      line = sprintf('  op_dat      args[%d] = {',nargs);
    elseif (l==2)
      word = 'idx'
      line = sprintf('  int         idxs[%d] = {',nargs);
    elseif (l==3)
      word = 'ptr'
      line = sprintf('  op_ptr      ptrs[%d] = {',nargs);
    elseif (l==4)
      word = 'dim'
      line = sprintf('  int         dims[%d] = {',nargs);
    elseif (l==5)
      word = 'typ'
      line = sprintf('  op_datatype typs[%d] = {',nargs);
    elseif (l==6)
      word = 'acc'
      line = sprintf('  op_access   accs[%d] = {',nargs);
    end

    for m = 1:nargs
      if (m<nargs)
        line = strcat(line,word,num2str(m-1),', ');
      else
        line = strcat(line,word,num2str(m-1),'};');
      end
    end
    file = strvcat(file,line);
  end

  line = sprintf('  int         inds[%d] = {',nargs);
  for m = 1:nargs
    if (m<nargs)
      line = strcat(line,num2str(inds(m)-1),', ');
    else
      line = strcat(line,num2str(inds(m)-1),'};');
    end
  end
  file = strvcat(file,line);


  file = strvcat(file,'  ',...
   '  if (OP_DIAGS>1) {              ',...
  ['    printf(" kernel routine with indirection: ' fn_name ' \n");'],...
   '  }                              ',' ',...
   '  // get plan                    ',' ',...
   '  op_plan *Plan = plan(name,set,nargs,args,idxs,ptrs,dims,typs,accs,ninds,inds);',' ',...
   '  // execute plan                ',' ',...
   '  int block_offset = 0;          ',' ');

  for m = 1:ninds
    if (indtyps(m)==OP_DOUBLE)
      line = sprintf('  double *ind_arg%d = arg%d.ddat_d;',m-1,invinds(m)-1);
    elseif (indtyps(m)==OP_FLOAT)
      line = sprintf('  float  *ind_arg%d = arg%d.fdat_d;',m-1,invinds(m)-1);
    elseif (indtyps(m)==OP_INT)
      line = sprintf('  int    *ind_arg%d = arg%d.idat_d;',m-1,invinds(m)-1);
    end
    file = strvcat(file,line);
  end

  file = strvcat(file,' ','  for (int col=0; col<(*Plan).ncolors; col++) { ',' ',...
                          '    int nblocks = (*Plan).ncolblk[col];         ',...
                          '    int nshared = (*Plan).nshared;              ',' ');

  line = ' ';
  for m = 1:nargs
    if (typs(m)==OP_DOUBLE)
      line = strcat(line,sprintf('arg%d.ddat_d,',m-1));
    elseif (typs(m)==OP_FLOAT)
      line = strcat(line,sprintf('arg%d.fdat_d,',m-1));
    elseif (typs(m)==OP_INT)
      line = strcat(line,sprintf('arg%d.idat_d,',m-1));
    end
  end
  file = strvcat(file,['    op_cuda_' fn_name '<<<nblocks,64,nshared>>>(']);

  for m = 1:ninds
    line = [ sprintf('        ind_arg%d,',          m-1) ...
             sprintf(' (*Plan).ind_ptrs[%d],', m-1) ...
             sprintf(' (*Plan).ind_sizes[%d],',m-1) ...
             sprintf(' (*Plan).ind_offs[%d],', m-1) ];
    file = strvcat(file,line);
  end

  for m = 1:nargs
    if (inds(m)==0)
      if (typs(m)==OP_DOUBLE)
        file = strvcat(file,sprintf('        arg%d.ddat_d,',m-1));
      elseif (typs(m)==OP_FLOAT)
        file = strvcat(file,sprintf('        arg%d.fdat_d,',m-1));
      elseif (typs(m)==OP_INT)
        file = strvcat(file,sprintf('        arg%d.idat_d,',m-1));
      end
    else
      file = strvcat(file,sprintf('        (*Plan).ptrs[%d],',m-1));
    end
  end

  file = strvcat(file, ... 
    '        block_offset,                                       ', ...
    '        (*Plan).blkmap,                                     ', ...
    '        (*Plan).offset,                                     ', ...
    '        (*Plan).nelems,                                     ', ...
    '        (*Plan).nthrcol,                                    ', ...
    '        (*Plan).thrcol);                                    ',' ', ...
    '    cutilSafeCall(cudaThreadSynchronize());                 ', ...
   ['    cutilCheckMsg("op_cuda_' fn_name ' execution failed\n");'],' ', ...
    '    block_offset += nblocks;                                ', ...
    '  }                                                         ', ...
    '}                                                           ');

  file


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% special version for kernels without indirect addressing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  else

  file = strvcat(file,strcat('__global__ void op_cuda_',fn_name,'('));
  for m = 1:nargs
    line = '                   TYP *ARG,';
    file = strvcat(file,rep(line,m));
  end
  file = strvcat(file,'                   int set_size) {',' ');

  for m = 1:nargs
    line = '  TYP ARG_l[DIM];';
    file = strvcat(file,rep(line,m));
  end

  file = strvcat(file,' ','  // process set elements',' ', ...
                      '  for (int n=threadIdx.x+blockIdx.x*blockDim.x;', ...
                      '       n<set_size; n+=blockDim.x*gridDim.x) {', ...
                      ' ','    // initialise local variables ',' ');

  for m = 1:nargs
    if (accs(m)==OP_READ | accs(m)==OP_RW)
      line = '    for (int d=0; d<DIM; d++)';
      file = strvcat(file,rep(line,m));
      line = '      ARG_l[d] = ARG[d+n*DIM];';
      file = strvcat(file,rep(line,m));

    elseif (accs(m)==OP_INC)
      line = '    for (int d=0; d<DIM; d++)';
      file = strvcat(file,rep(line,m));
      line = '      ARG_l[d] = 0;';
      file = strvcat(file,rep(line,m));
    end
  end

  file = strvcat(file,' ','    // user-supplied kernel call',' ');

  line = regexprep('    FN(','FN',fn_name);
  for m = 1:nargs
    line = strcat(line,rep('ARG_l,',m));
  end
  line = strcat(line(1:end-1),');');
  file = strvcat(file,line,' ');

  file = strvcat(file,'    // store local variables ',' ');

  for m = 1:nargs
    if (accs(m)==OP_WRITE | accs(m)==OP_RW)
      line = '    for (int d=0; d<DIM; d++)';
      file = strvcat(file,rep(line,m));
      line = '      ARG[d+n*DIM]  = ARG_l[d];';
      file = strvcat(file,rep(line,m));

    elseif (accs(m)==OP_INC)
      line = '    for (int d=0; d<DIM; d++)';
      file = strvcat(file,rep(line,m));
      line = '      ARG[d+n*DIM] += ARG_l[d];';
      file = strvcat(file,rep(line,m));
    end
  end

  file = strvcat(file,'  }','}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% add stub function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  file = strvcat(file,' ',' ',...
                      '// host stub function            ',' ',...
	         'extern "C"                       ',...
                ['void op_par_loop_' fn_name '(char const * name, op_set set,']);

  for m = 1:nargs
    line = '  op_dat ARG, int ARGidx, op_ptr ARGptr, int ARGdim,';
    file = strvcat(file,rep(line,m));
    if (m==nargs)
      line = '       op_datatype ARGtyp,           op_access ARGacc){';
    else
      line = '       op_datatype ARGtyp,           op_access ARGacc,';
    end
    file = strvcat(file,rep(line,m));
  end

  file = strvcat(file,'  ','  if (OP_DIAGS>1) {              ',...
                     ['    printf(" kernel routine w/o indirection:  ' fn_name ' \n");'],...
                      '  }                              ',' ',...
                      '  // execute plan                ',' ');
  line = ' ';
  for m = 1:nargs
    if (typs(m)==OP_DOUBLE)
      line = strcat(line,rep('ARG.ddat_d,',m));
    elseif (typs(m)==OP_FLOAT)
      line = strcat(line,rep('ARG.fdat_d,',m));
    elseif (typs(m)==OP_INT)
      line = strcat(line,rep('ARG.idat_d,',m));
    end
  end
  file = strvcat(file,['  op_cuda_' fn_name '<<<100,64>>>(' line 'set.size);']);
  file = strvcat(file,' ','  cutilSafeCall(cudaThreadSynchronize());');
  file = strvcat(file,['  cutilCheckMsg("op_cuda_', fn_name ' execution failed\n");']);
  file = strvcat(file,'}                                ');

  %file

  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% output to a CUDA file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  fid = fopen(strcat(fn_name,'_kernel.cu'),'wt');
  fprintf(fid,'// \n// auto-generated by op2.m on %s \n//\n\n',datestr(now));
  for n=1:size(file,1)
    fprintf(fid,'%s\n',file(n,:));
  end
  fclose(fid);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% finally, output new main file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


src_file = fileread([filename '.cpp']);

code_segment = ' ';

for k = 1:nkernels
  nargs = kernel_args(k);
  old = [ 'op_par_loop_' num2str(nargs) '(' kernel{k} ','];
  new = [ 'op_par_loop_' kernel{k} '('];
  src_file = regexprep(src_file,old,new);

  code_segment = strvcat(code_segment,...
   'extern "C"                                             ' ,...
  ['void op_par_loop_' kernel{k} '(char const *, op_set,   ']);

  for n = 1:nargs-1
    code_segment = strvcat(code_segment,...  
       '  op_dat, int, op_ptr, int, op_datatype, op_access,    ');
  end
    code_segment = strvcat(code_segment,...  
       '  op_dat, int, op_ptr, int, op_datatype, op_access);   ',' ');
end


code_segment = strvcat( ...
'#include "op_datatypes.h"                                                             ',...
'                                                                                      ',... 
'extern "C"                                                                            ',...
'void op_init(int, char **);                                                           ',...
'                                                                                      ',...
'extern "C"                                                                            ',...
'void op_decl_set(int, op_set &, char const *);                                        ',...
'                                                                                      ',... 
'extern "C"                                                                            ',... 
'void op_decl_ptr(op_set, op_set, int, int *, op_ptr &, char const *);                 ',...
'                                                                                      ',... 
'extern "C"                                                                            ',... 
'void op_decl_ddat(op_set, int, op_datatype, double *, op_dat &, char const *);        ',...
'                                                                                      ',... 
'extern "C"                                                                            ',... 
'void op_decl_fdat(op_set, int, op_datatype, float *, op_dat &, char const *);         ',...
'                                                                                      ',... 
'extern "C"                                                                            ',...  
'void op_decl_idat(op_set, int, op_datatype, int *, op_dat &, char const *);           ',...
'                                                                                      ',... 
'void op_decl_dat(op_set s,int dm,op_datatype t,double *dat,op_dat &dt,char const *nm){',...
'    op_decl_ddat(       s,    dm,            t,        dat,        dt,            nm);',...
'}                                                                                     ',... 
'                                                                                      ',... 
'void op_decl_dat(op_set s,int dim,op_datatype t,float *dat,op_dat &dt,char const *nm){',...
'    op_decl_fdat(       s,    dim,            t,       dat,        dt,            nm);',...
'}                                                                                     ',... 
'                                                                                      ',... 
'void op_decl_dat(op_set s,int dim,op_datatype t,int *dat,op_dat &data,char const *nm){',...
'    op_decl_idat(       s,    dim,            t,     dat,        data,            nm);',...
'}                                                                                     ',... 
'                                                                                      ',... 
'extern "C"                                                                            ',... 
'void op_decl_dconst(int, op_datatype, double *, char const *);                        ',...
'                                                                                      ',... 
'extern "C"                                                                            ',... 
'void op_decl_fconst(int, op_datatype, float *, char const *);                         ',...
'                                                                                      ',... 
'extern "C"                                                                            ',...  
'void op_decl_iconst(int, op_datatype, int *, char const *);                           ',...
'                                                                                      ',... 
'void op_decl_const(int dim,op_datatype type,double *dat,char const *nm){              ',...
'    op_decl_dconst(    dim,            type,        dat,            nm);              ',...
'}                                                                                     ',... 
'                                                                                      ',... 
'void op_decl_const(int dim,op_datatype type,float *dat,char const *nm){               ',...
'    op_decl_fconst(    dim,            type,       dat,            nm);               ',...
'}                                                                                     ',... 
'                                                                                      ',... 
'void op_decl_const(int dim,op_datatype type,int *dat,char const *nm){                 ',...
'    op_decl_iconst(    dim,            type,     dat,            nm);                 ',...
'}                                                                                     ',... 
'                                                                                      ',...  
'extern "C"                                                                            ',... 
'void op_fetch_data(op_dat);                                                           ',... 
'                                                                                      ',... 
'extern "C"                                                                            ',... 
'void op_diagnostic_output();                                                          ',... 
'                                                                                      ',... 
'//                                                                                    ',... 
'// op_par_loop declarations                                                           ',... 
'//                                                                                    ',... 
code_segment);


loc = strfind(src_file,'#include "op_seq.h"');

fid = fopen(strcat(filename,'_op.cpp'),'wt');
fprintf(fid,'// \n// auto-generated by op2.m on %s \n//\n\n',datestr(now));
fprintf(fid,'%s',src_file(1:loc-1));

for n=1:size(code_segment,1)
  fprintf(fid,'%s\n',code_segment(n,:));
end

fprintf(fid,'%s',src_file(loc+20:end));

fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ... and master kernel file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


file = strvcat('// header files          ',...
               '                         ',...
               '#include <stdlib.h>      ',...
               '#include <stdio.h>       ',...
               '#include <string.h>      ',...
               '#include <math.h>        ',...
               '#include <cutil_inline.h>',...
               '#include "op_datatypes.h"',...
               '                         ',...
               '// global constants      ',' ');

src_file = fileread([filename '.cpp']);
src_file = regexprep(src_file,'\s','');

while (~isempty(strfind(src_file,'op_decl_const(')))
  loc  = min(strfind(src_file,'op_decl_const('));
  src_file = src_file(loc+14:end);
  [src_args, src_file] = strtok(src_file,')');

  loc = [0 strfind(src_args,',') length(src_args)+1];
  na  = length(loc)-1;

  if( na ~= 4)
    disp(sprintf('wrong number of arguments in op_decl_const'));
    error('aborting')
  end

  for n = 1:na
    C{n} = src_args(loc(n)+1:loc(n+1)-1);
  end

  switch C{2}
    case 'OP_FLOAT'
      type = 'float';
    case 'OP_DOUBLE'
      type = 'double';
    case 'OP_INT'
      type = 'int';
    otherwise
      disp('unknown data type')
  end

  dim = str2num(C{1});
  if (dim==1)
    file = strvcat(file,[ '__constant__ ' type ' ' C{4}(2:end-1) ';' ]);
  else
    file = strvcat(file,[ '__constant__ ' type ' ' C{4}(2:end-1) '[' C{1} '];' ]);
  end

end

  file = strvcat(file,' ','// user kernel files',' ');

for k = 1:nkernels
  file = strvcat(file, ['#include "' kernel{k} '_kernel.cu"']);
end



fid = fopen(strcat(filename,'_kernels.cu'),'wt');
fprintf(fid,'// \n// auto-generated by op2.m on %s \n//\n\n',datestr(now));

for n=1:size(file,1)
  fprintf(fid,'%s\n',file(n,:));
end
fclose(fid);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% a little function to replace keywords
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function line = rep(line,m)

global dims idxs typs indtyps inddims
global OP_FLOAT OP_DOUBLE OP_INT

typ{OP_FLOAT}  = 'float ';
typ{OP_DOUBLE} = 'double';
typ{OP_INT}    = 'int   ';

line = regexprep(line,'INDDIM',num2str(inddims(m)));
line = regexprep(line,'INDARG',sprintf('ind_arg%d',m-1));
line = regexprep(line,'INDTYP',typ{indtyps(m)});

line = regexprep(line,'DIM',num2str(dims(m)));
line = regexprep(line,'ARG',sprintf('arg%d',m-1));
line = regexprep(line,'TYP',typ{typs(m)});
line = regexprep(line,'IDX',num2str(idxs(m)));
