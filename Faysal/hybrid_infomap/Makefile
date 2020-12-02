# Various flags
CXX  = CC
LINK = $(CXX)
#CXXFLAGS = -I -Wall -g 
CXXFLAGS = -g -Wall -O0 -fopenmp  #-I #-Wall -O3 -funroll-loops -pipe 
#CXXFLAGS = -g -Wall -fopenmp #-I #-Wall -O3 -funroll-loops -pipe 
LFLAGS =  -g -fopenmp -Wall -Werror -O3

TARGET  = ompInfomap

HEADER  = Node.h Module.h FileIO.h timing.h
FILES = OmpRelaxmap.cpp Node.cpp Module.cpp FileIO.cpp timing.cpp

OBJECTS = $(FILES:.cpp=.o)

$(TARGET): ${OBJECTS}
	$(LINK) $(LFLAGS) -std=c++11 $^ -o $@ -lmetis 

all: $(TARGET)

clean:
	rm -f $(OBJECTS)

distclean:
	rm -f $(OBJECTS) $(TARGET)

# Compile and dependency
$(OBJECTS): $(HEADER) Makefile




