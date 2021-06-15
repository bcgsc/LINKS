GCC := g++

OUTPUT := LINKS 
SOURCES := $(wildcard *.hpp) $(wildcard *.cpp)  
CCFLAGS := -std=c++11 -pthread -I../../../btllib/include/

all: $(OUTPUT)

$(OUTPUT):
	$(GCC) -o $(OUTPUT) $(CCFLAGS) $(SOURCES)

clean:
	rm $(OUTPUT)

.PHONY: all
