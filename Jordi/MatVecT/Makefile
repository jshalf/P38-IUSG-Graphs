SRCDIR = src
CPPFILES = \
	$(SRCDIR)/Misc.cpp \
	$(SRCDIR)/MatVec.cpp \
	$(SRCDIR)/Matrix.cpp \
	$(SRCDIR)/Main.cpp \

COMPILE = g++ -fopenmp -lm -O0 -std=c++11

all: clean Main

Main: $(SRCDIR)/Main.cpp
	$(COMPILE) $(CPPFILES) -o Main

clean:
	rm -f Main
