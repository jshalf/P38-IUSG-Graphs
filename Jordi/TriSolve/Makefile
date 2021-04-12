SRCDIR = ../src
TRISOLVE_SRCDIR = src

CPPFILES = \
	$(SRCDIR)/Misc.cpp \
	$(SRCDIR)/Matrix.cpp \
        $(SRCDIR)/MsgQ.cpp \
        $(TRISOLVE_SRCDIR)/TriSolve.cpp \
	$(TRISOLVE_SRCDIR)/Main.cpp \

COMPILE = g++ -fopenmp -lm -O3 -std=c++11 -g

LIBS =
INCLUDE =

all: clean Main

Main: $(TRISOLVE_SRCDIR)/Main.cpp
	$(COMPILE) $(CPPFILES) $(INCLUDE) -o Main $(LIBS)

clean:
	rm -f Main
