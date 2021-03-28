# Various flags
CXX  = CC
LINK = $(CXX)
#CXXFLAGS = -I -Wall -g 
CXXFLAGS = -g -Wall -O0 -fopenmp  #-I #-Wall -O3 -funroll-loops -pipe 
#CXXFLAGS = -g -Wall -fopenmp #-I #-Wall -O3 -funroll-loops -pipe 
LFLAGS =  -g -fopenmp -Wall -Werror -O3

VTUNE_LIB_DIR =/usr/common/software/intel/parallel_studio_xe_2020_cluster_edition/vtune_profiler/lib64/libittnotify.a 
INCADD =-I/opt/intel/vtune_profiler_2020.2.0.610396/include


TARGET  = ompInfomap

HEADER  = Node.h Module.h FileIO.h timing.h global.h
FILES = OmpRelaxmap.cpp Node.cpp Module.cpp FileIO.cpp timing.cpp global.cpp

OBJECTS = $(FILES:.cpp=.o)

$(TARGET): ${OBJECTS}
	$(LINK) $(LFLAGS) -std=c++11 $^ $(INCADD) $(VTUNE_LIB_DIR) -o $@ -lmetis 

all: $(TARGET)

clean:
	rm -f $(OBJECTS)

distclean:
	rm -f $(OBJECTS) $(TARGET)

# Compile and dependency
$(OBJECTS): $(HEADER) Makefile




