NVCC=$(CUDA_HOME)/bin/nvcc
NVCCFLAGS = -O3 -G --shared -Xcompiler -fPIC $(R_INCL) -arch=sm_30 --relocatable-device-code true -I$(CUDA_HOME)/include
NVCCLINKS = -shared -lcuda -lcublas -lcudadevrt -Wno-deprecated-gpu-targets
#OBJECTS = dev-fnc.o kernel.o
#PKG_LIBS = $(OBJECTS) -L$(CUDA_HOME)/lib64 -lcuda -lcublas -lcudadevrt
cuda_obj = dev-fnc.o kernel.o

all: cuda

cuda:
	$(NVCC) $(NVCCFLAGS) -c dev-fnc.cu -o dev-fnc.o
	$(NVCC) $(NVCCFLAGS) -c kernel.cu -o kernel.o
	$(NVCC) $(NVCCLINKS) $(cuda_obj) -o uroot.so
