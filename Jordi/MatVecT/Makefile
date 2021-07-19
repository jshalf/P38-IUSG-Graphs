SRCDIR = ../src
MATVECT_SRCDIR = src

CPPFILES = \
	$(SRCDIR)/Misc.cpp \
	$(SRCDIR)/Matrix.cpp \
	$(SRCDIR)/MsgQ.cpp \
        $(MATVECT_SRCDIR)/MatVecT.cpp \
	$(MATVECT_SRCDIR)/Main.cpp \

COMPILE = g++ -fopenmp -lm -O3 -std=c++11

all: clean Main

Main: $(MATVECT_SRCDIR)/Main.cpp
	$(COMPILE) $(CPPFILES) -o Main

clean:
	rm -f Main
