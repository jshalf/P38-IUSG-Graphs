SRCDIR = ../src
ASYNC_JACOBI_SRCDIR = src

CPPFILES = \
	$(SRCDIR)/Misc.cpp \
	$(SRCDIR)/Matrix.cpp \
        $(SRCDIR)/MsgQ.cpp \
	$(ASYNC_JACOBI_SRCDIR)/Jacobi_CSR.cpp \
	$(ASYNC_JACOBI_SRCDIR)/Main.cpp \

COMPILE = g++ -fopenmp -lm -O3 -std=c++11

all: clean Main

Main: $(ASYNC_JACOBI_SRCDIR)/Main.cpp
	$(COMPILE) $(CPPFILES) -o Main

clean:
	rm -f Main
