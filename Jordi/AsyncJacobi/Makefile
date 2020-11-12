SRCDIR = src
CPPFILES = \
	$(SRCDIR)/Misc.cpp \
	$(SRCDIR)/Matrix.cpp \
	$(SRCDIR)/Jacobi.cpp \
	$(SRCDIR)/Main.cpp \

COMPILE = g++ -fopenmp -lm -O3 -std=c++11

all: clean Main

Main: $(SRCDIR)/Main.cpp
	$(COMPILE) $(CPPFILES) -o Main

clean:
	rm -f Main
