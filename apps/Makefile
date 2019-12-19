# Change prefix and/or ROOTDIR, if required.
PREFIX = /localdisk
ROOTDIR=/wrk3

# Other
CPPUTIL = $(ROOTDIR)/simploce/pt-cgmd/cpputil
PARTICLES = $(ROOTDIR)/simploce/pt-cgmd/particles
SIMULATION = $(ROOTDIR)/simploce/pt-cgmd/simulation

# Relative to build directory.
VPATH = ../src

# All sources.
APPS = s-cgmd.cpp s-gr.cpp s-dipole-moment.cpp s-msd.cpp

# All executables.
AEXE = $(APPS:.cpp=)

# For g++ >= 5.4
CC = g++
LT = /usr/bin/libtool
OPT = -ggdb
STD = c++14
INCLPATH = -I$(CPPUTIL)/include -I$(PARTICLES)/include -I$(SIMULATION)/include

# Libraries.
OLIBS = -lcpputil -lparticles -lsimulation

# Compiler and linker flags
CFLAGS = $(INCLPATH) $(OPT) -Wall -std=$(STD) -pthread -D _DEBUG
LDFLAGS =  $(INCLPATH) $(OPT) -Wall -std=$(STD) -L. -L$(PREFIX)/lib -pthread
LIBS = -lm -lpthread $(OLIBS) -lboost_program_options -lboost_iostreams
# -lbz2 -lz
#-lz for zlib

# All object files.
%.o : %.cpp
	$(LT) --mode=compile $(CC) -c $(CFLAGS) $< -o $@

# All executables.
% : %.cpp
	$(LT) --mode=link $(CC) $(LDFLAGS) $(LIBS) $< -o $@ -rpath $(PREFIX)/bin

apps: $(AEXE)
	echo "Done: apps"

all: apps
	echo "Done: all"

install: all
	$(LT) --mode=install install $(AEXE) $(PREFIX)/bin

clean-build :
	$(LT) --mode=clean \rm -r -f ../build/.libs/*
	$(LT) --mode=clean \rm -r -f ../build/.libs/*.*
	\rm -r -f ../build/.libs
	$(LT) --mode=clean \rm -r -f ../build/*
	echo "Done: clean-build"

clean-apps:
	$(LT) --mode=clean \rm -r -f $(AEXE)
	echo "Done: clean-apps"

clean: clean-apps clean-build
	echo "Done: clean"
