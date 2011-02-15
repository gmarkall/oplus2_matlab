################################################################################
#
# Build script for project
#
################################################################################

# Executable
EXECUTABLE	:= jac

# CUDA source files (compiled with nvcc)
CUFILES		:= jac_kernels.cu  op_lib.cu

# CUDA dependency files
CU_DEPS		:= op_datatypes.h \
		   res_kernel.cu \
		   update_kernel.cu 

# C/C++ source files (compiled with gcc / c++)
CCFILES		:= jac_op.cpp

# C/C++ dependency files
CC_DEPS		:= op_datatypes.h


################################################################################
# Rules and targets

ROOTDIR = ../NVIDIA_CUDA_SDK_2.2/common
BINDIR  = ./bin
include $(ROOTDIR)/common.mk

NVCCFLAGS += -Xptxas -v

NVCCFLAGS += -arch sm_13
