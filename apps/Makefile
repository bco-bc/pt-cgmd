PREFIX = /home/ajuffer

# Relative to apps/build folder.
CPPUTIL = $(PREFIX)/simploce/pt-cgmd/cpputil
PARTICLES = $(PREFIX)/simploce/pt-cgmd/particles
SIMULATION = $(PREFIX)/simploce/pt-cgmd/simulation

VPATH = ../src

APPS = s-cgmd.cpp s-gr.cpp

AEXE = $(APPS:.cpp=)

# For g++ >= 5.4
CC = g++
LT = libtool
OPT = -ggdb
STD = c++14
INCLPATH = -I$(CPPUTIL)/include -I$(PARTICLES)/include -I$(SIMULATION)/include

OLIBS = -lcpputil -lparticles -lsimulation

CFLAGS = $(INCLPATH) $(OPT) -Wall -std=$(STD) -pthread
LDFLAGS =  $(INCLPATH) $(OPT) -Wall -std=$(STD) -L. -L$(PREFIX)/lib -pthread
LIBS = -lm -lpthread $(OLIBS) -lboost_program_options -lboost_iostreams
# -lbz2 -lz
#-lz for zlib

# All object files. Default.
%.o : %.cpp
	$(LT) --mode=compile $(CC) -c $(CFLAGS) $< -o $@

# All executables.
% : %.cpp $(INCLUDE)
	$(LT) --mode=link $(CC) $(LDFLAGS) $(LIBS) $< -o $@ -rpath $(PREFIX)/bin

apps: $(AEXE)
	echo "Done: apps"

all: apps
	echo "Done: all"

install: all
	$(LT) --mode=install install $(AEXE) $(PREFIX)/bin

doc: $(INCLUDE) $(TEMPLATES)
	cd ..;pwd;/usr/bin/doxygen Doxyfile;

clean-build :
	$(LT) --mode=clean \rm -r -f ../build/.libs/*
	$(LT) --mode=clean \rm -r -f ../build/.libs/*.*
	\rm -r -f ../build/.libs
	$(LT) --mode=clean \rm -r -f ../build/*

clean-apps:
	$(LT) --mode=clean \rm -r -f $(AEXE)

clean: clean-apps clean-build
	echo "Done: clean"
