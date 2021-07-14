CXX = icpc
CXXFLAGS = -std=c++11 -ipo -O3 -static -fp-model fast=2 -xHost
LDFLAGS = -std=c++11
MKLFLAGS = -mkl=parallel
NVCC = nvcc
NVFLAGS= -std=c++11 -ccbin=icpc -use_fast_math -arch=sm_60
NVLIBS= -lcudart
NVLIBDIRS=-L/usr/local/cuda/lib64
NVINCDIRS=-I/usr/local/cuda/include

OBJS = obj/main.o obj/io.o obj/energy.o obj/xyz.o obj/moves.o obj/gpu.o
HEADERS = conf.h io.h energy.h xyz.h moves.h gpu.h

mc.exe:	$(OBJS) $(HEADERS)
	$(CXX) $(LDFLAGS) $(MKLFLAGS) ${NVLIBS} ${NVLIBDIRS} ${NVINCDIRS} -o ../$@ $^

obj/main.o: main.c $(HEADERS)
	$(CXX) -c $(CXXFLAGS) -o $@ main.c

obj/%.o: %.c $(HEADERS)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

obj/gpu.o: gpu.cu $(HEADERS)
	${NVCC} ${NVFLAGS} ${NVLIBS} -Xcompiler "${MKLFLAGS}" -c -o $@ $<

clean:
	rm -fr ../mc.exe $(OBJS)
