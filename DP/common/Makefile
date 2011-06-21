
MPI_HOME	:= /home/gihan/openmpi

CC		:=$(MPI_HOME)/bin/mpicc 
CXX		:=$(MPI_HOME)/bin/mpiCC
CCFLAGS		:= -O2  -vec-report -xSSE2,SSE3,SSE4.1,SSE4.2  -g #-DCOMM_PERF #-DDEBUG


#
# ParMetis for partitioning
#
PARMETIS_HOME	:= /home/gihan/ParMetis-3.1.1
PARMETISLIB 	:= $(PARMETIS_HOME)/libparmetis.a $(PARMETIS_HOME)/libmetis.a
INC		:= $(INC) -I$(PARMETIS_HOME)  


#
# make all by default
#

all:		op_lib_core.o op_lib.o op_lib.a

#
# core OP2 library written in C99 for interoperability with FORTRAN 9x
#

op_lib_core.o:	op_lib_core.c op_lib_core.h op_mpi_core.c op_mpi_core.h \
		$(PARMETIS_LIB) op_mpi_part_core.c op_mpi_part_core.h mpi_util.c
		$(CC) $(CCFLAGS) $(INC) -std=c99 -c op_lib_core.c 
		

#
# C++ library used by sequential and OpenMP executables
#

op_lib.o:	op_lib.cpp op_lib.h op_lib_core.h
		$(CXX) $(CCFLAGS) -c op_lib.cpp


op_lib.a:	op_lib_core.o op_lib.o
		rm -f op_lib.a
		ar rcs op_lib.a op_lib_core.o op_lib.o 

#
# clean up
#

clean:
		rm -f op_lib_core.o op_lib.o op_lib.a
