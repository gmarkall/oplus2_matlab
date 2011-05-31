
#
# make all by default
#

all:		op_lib_core.o op_lib.o op_lib.a

#
# core OP2 library written in C99 for interoperability with FORTRAN 9x
#

op_lib_core.o:	op_lib_core.c op_lib_core.h
		gcc -g -std=c99 -c op_lib_core.c

#
# C++ library used by sequential and OpenMP executables
#

op_lib.o:	op_lib.cpp op_lib.h op_lib_core.h
		g++ -g -c op_lib.cpp


op_lib.a:	op_lib_core.o op_lib.o
		ar rcs op_lib.a op_lib_core.o op_lib.o

#
# clean up
#

clean:
		rm -f op_lib_core.o op_lib.o op_lib.a
