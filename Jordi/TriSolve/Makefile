SHELL = bash

CXX = g++
PP_FLAGS = -O0 -std=c++14
CG_FLAGS =
LD_FLAGS =
LIBS =

SRCDIR = ../src
TRISOLVE_SRCDIR = src

CPPFILES = \
	$(SRCDIR)/Misc.cpp \
	$(SRCDIR)/Matrix.cpp \
	$(SRCDIR)/LevelSets.cpp \
        $(SRCDIR)/Parallel.cpp \
        $(TRISOLVE_SRCDIR)/TriSolve.cpp \
	$(TRISOLVE_SRCDIR)/Main.cpp \

INCLUDE = -I../src -Isrc
LIBS = -fopenmp -lm

#Configure parallel backbend can be one of the following
# {omp, pthreads, deva}
parallel_backend ?= omp

ifeq ($(par_backend),deva)
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

  CXX = $(deva_cxx)
  PP_FLAGS += -DUSE_DEVA $(deva_ppflags)
  CG_FLAGS += $(deva_cgflags)
  LD_FLAGS += $(deva_ldflags)
  #this is a slight hack as it appends compile devastator objects to a list of source files
  CPPFILES += $(deva_objs)
  LIBS += $(deva_libflags)
else ifeq ($(parallel_backend),pthreads)
  PP_FLAGS += -DUSE_PTHREADS
else ifeq ($(parallel_backend),omp)
  #no-op
else
  $(error Invalid parallel backend of type: $(parallel_backend). Must be "omp", "deva", or "pthreads")
endif

all: Main

Main: $(CPPFILES)
	$(CXX) $(PP_FLAGS) $(INCLUDE) $(CG_FLAGS) $(LD_FLAGS) -o "$@" $(CPPFILES) $(LIBS)

clean:
	rm -f Main

.PHONY: all clean
