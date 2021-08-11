GCC := /gsc/software/linux-x86_64-centos7/gcc-7.2.0/bin/g++

OUTPUT := LINKS 
SOURCES := $(wildcard *.hpp) $(wildcard *.cpp)  
CCFLAGS := -std=c++17 -pthread  -I../../../btllib/include/ -I/gsc/software/linux-x86_64-centos7/boost_1_55_0/ -L/projects/btl_scratch/tgoktas/miniconda3/lib/ 
#CCFLAGS := -std=c++17 -lunwind -pthread  -I../../../btllib/include/ -L/projects/btl_scratch/tgoktas/miniconda3/lib/ -lprofiler

all: $(OUTPUT)

$(OUTPUT):
	$(GCC) -o $(OUTPUT) $(CCFLAGS) $(SOURCES)

clean:
	rm $(OUTPUT)

.PHONY: all
