#required an empty 'uroot.so' object even if CUDA not available 
#because NAMESPACE requires useDynLib(uroot);
#do not use an empty 'uroot.cpp' file, otherwise default Makefile is run, 
#which could be used to do the linking using the options defined there but 
#gave error "undefined symbol: cudaFree", therefore linking is done using 'nvcc' 
#in COMPILE_CUDA_TRUE.mk
#
#g++ -shared -L/usr/lib/R/lib -lR -o uroot.so

all: cuda

cuda:
	$(CXX) -shared $(LDFLAGS) -lR -o uroot.so
