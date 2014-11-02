#
# set paths for header files and libraries
#

OP2		:= ../common
CUDA		:= $(CUDA_INSTALL_PATH)
INC		:= -I$(CUDA)/include -I$(OP2) -I.
LIB		:= -L$(CUDA)/lib64

#
# set flags for compilation and linking
#

CCFLAGS		:= -g -fPIC -DUNIX -Wall -Wextra
NVCCFLAGS	:= -g -G -arch=sm_20 -Xptxas=-v -use_fast_math
LIBS		:= $(LIB) -lcudart
OMP		:= -fopenmp -lgomp

#
# master to make all versions
#

all:		test_seq

#
# simple sequential version
#

test_seq:	test.cpp $(OP2)/op_seq.cpp $(OP2)/op_datatypes.h
		g++ $(CCFLAGS) $(INC) -o test_seq test.cpp $(OP2)/op_seq.cpp

clean:
		rm -f test_seq
