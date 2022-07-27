SHELL = bash

# assumes devastator is installed in sibling directory

app_dir := $(shell pwd)
deva_dir := $(shell cd ..; pwd)/devastator

# canonicalize devastator configuration parameters
world ?= threads
ifeq ($(world),threads)
  ranks ?= 1
else ifeq ($(world),gasnet)
  procs ?= 32
  workers ?= 31
else
  $(error unknown world: $(world))
endif
debug ?= 0
ifneq ($(debug),0)
  optlev ?= 0
  syms ?= 1
  opnew ?= libc
  asan ?= 1
else
  optlev ?= 3
  syms ?= 0
  opnew ?= deva
  asan ?= 0
endif
dummy ?= 0

# populate this with configuration variables for devastator, this logic uses
# the same names devastator expects, but doesn't have to, e.g. the app could
# expect "foo=..." from the user and give that to devastator as "bar=$(foo)"
ifeq ($(world), threads)
  deva_conf := world=threads ranks=$(ranks)
else ifeq ($(world), gasnet)
  deva_conf := world=gasnet procs=$(procs) workers=$(workers)
endif
deva_conf := $(deva_conf) \
             debug=$(debug) \
             optlev=$(optlev) \
             syms=$(syms) \
             opnew=$(opnew) \
             asan=$(asan) \
             dummy=$(dummy)

# function for invoking "brutal $(deva_conf) makevar $(1) $(2)"
deva_var = $(shell cd $(deva_dir); . sourceme; cd $(app_dir); DEVA=$(deva_dir) brutal $(deva_conf) makevar $(1) $(2))

deva_inclhpp_dir = src

# extract build flags from devastator given that deva_inclhpp_dir/deva_includes.hpp
# includes all the functionality this app needs from devastator
deva_cxx      := $(call deva_var, CXX, $(deva_inclhpp_dir)/deva_includes.hpp)
deva_ppflags  := $(call deva_var, PPFLAGS, $(deva_inclhpp_dir)/deva_includes.hpp)
deva_cgflags  := $(call deva_var, CGFLAGS, $(deva_inclhpp_dir)/deva_includes.hpp)
deva_ldflags  := $(call deva_var, LDFLAGS, $(deva_inclhpp_dir)/deva_includes.hpp)
deva_objs     := $(call deva_var, OBJS, $(deva_inclhpp_dir)/deva_includes.hpp)
deva_libflags := $(call deva_var, LIBFLAGS, $(deva_inclhpp_dir)/deva_includes.hpp)

SRCDIR = ../src
TRISOLVE_SRCDIR = src

CPPFILES = \
	$(SRCDIR)/Misc.cpp \
	$(SRCDIR)/Matrix.cpp \
	$(SRCDIR)/LevelSets.cpp \
        $(SRCDIR)/Parallel.cpp \
        $(TRISOLVE_SRCDIR)/TriSolve.cpp \
	$(TRISOLVE_SRCDIR)/Main.cpp \

PP_FLAGS = -O0 -std=c++14 #-DUSE_DEVA #-DUSE_PTHREADS
INCLUDE = -I../src -Isrc
LIBS = -fopenmp -lm 

all: Main

Main: $(CPPFILES) $(deva_objs)
	$(deva_cxx) $(PP_FLAGS) $(INCLUDE) $(deva_ppflags) $(deva_cgflags) $(deva_ldflags) -o "$@" $(deva_objs) $(CPPFILES) $(LIBS) $(deva_libflags)

clean:
	rm -r $(deva_dir)/.brutal
	rm -f Main

.PHONY: all clean
