SRCDIR = ../src
TRISOLVE_SRCDIR = src

CPPFILES = \
	$(SRCDIR)/Misc.cpp \
	$(SRCDIR)/Matrix.cpp \
	$(TRISOLVE_SRCDIR)/TriSolve_COO.cpp \
        $(TRISOLVE_SRCDIR)/TriSolve_CSR.cpp \
	$(TRISOLVE_SRCDIR)/Main.cpp \

COMPILE = g++ -fopenmp -lm -O0 -std=c++11 -g

LIBS = $(SUPERLU_DIR)/lib/libsuperlu_mt_OPENMP.a $(SUPERLU_DIR)/lib/libblas_OPENMP.a
INCLUDE = -I$(SUPERLU_DIR)/SRC/

all: clean Main

SUPERLU_DIR = ../SuperLU_MT_3.1

Main: $(TRISOLVE_SRCDIR)/Main.cpp
	$(COMPILE) $(CPPFILES) $(INCLUDE) -o Main $(LIBS)

clean:
	rm -f Main
