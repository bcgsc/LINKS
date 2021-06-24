GCC := /gsc/software/linux-x86_64-centos7/gcc-7.2.0/bin/g++

OUTPUT := LINKS 
SOURCES := $(wildcard *.hpp) $(wildcard *.cpp)  
CCFLAGS := -std=c++17 -pthread -I../../../btllib/include/

all: $(OUTPUT)

$(OUTPUT):
	$(GCC) -o $(OUTPUT) $(CCFLAGS) $(SOURCES)

clean:
	rm $(OUTPUT)

.PHONY: all
