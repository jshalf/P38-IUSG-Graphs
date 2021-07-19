SRCDIR = ../src
ATOMIC_SRCDIR = src

CPPFILES = \
	$(SRCDIR)/Misc.cpp \
        $(SRCDIR)/MsgQ.cpp \
	$(ATOMIC_SRCDIR)/Main.cpp \

COMPILE = g++ -fopenmp -lm -O3 -std=c++11

all: clean Main

Main: $(ATOMIC_SRCDIR)/Main.cpp
	$(COMPILE) $(CPPFILES) -D WINDOWS -o Main

clean:
	rm -f Main
