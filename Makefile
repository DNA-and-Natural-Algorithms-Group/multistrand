#
# Copyright (c) 2007-2010 Caltech. All rights reserved.
#   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)   [maintainer]
#  	             Chris Berlind (cberlind@dna.caltech.edu)
#

# A useful line for debugging dependencies:
#    @echo Rule expansion: $@ is target, [$?] deps newer, all: [$^].

vpath %.h include
INCLUDES = energymodel.h loop.h move.h optionlists.h options.h python_options.h scomplex.h scomplexlist.h ssystem.h strandordering.h

SOURCES_LOOP = loop.cc move.cc
SOURCES_ENERGYMODEL = energymodel.cc nupackenergymodel.cc viennaenergymodel.cc
SOURCES_STATE = scomplex.cc scomplexlist.cc strandordering.cc
SOURCES_SYSTEM = ssystem.cc
SOURCES_OPTIONS = optionlists.cc python_options.cc
SOURCES_TESTING = embedding_test.cc embedding_test2.cc testingmain.cc

SOURCES = $(SOURCES_LOOP) $(SOURCES_ENERGYMODEL) $(SOURCES_STATE) $(SOURCES_SYSTEM) $(SOURCES_OPTIONS) $(SOURCES_TESTING)

# separate the objects to make it a bit easier on ourselves when compiling targets
#  for the python library interface.

VPATH=loop state system energymodel include interface

MAIN_OBJECT = testingmain.o
OBJECTS = loop.o scomplex.o energymodel.o viennaenergymodel.o nupackenergymodel.o move.o ssystem.o scomplexlist.o strandordering.o python_options.o optionlists.o
INTERFACE_OBJECTS := boost_interface.o $(OBJECTS)

OBJPATH=obj
OBJECTS := $(OBJECTS:%=$(OBJPATH)/%)
MAIN_OBJECT := $(MAIN_OBJECT:%=$(OBJPATH)/%)
INTERFACE_OBJECTS := $(INTERFACE_OBJECTS:%=$(OBJPATH)/interface/%)


MULTISTRAND_INCLUDES := -I $(realpath ./include)

PYTHON_INCLUDES 	  = -I /usr/include/python2.6
BOOST_INCLUDES 	      = -I /opt/local/include
# CHANGE THESE TWO TO REFLECT YOUR PATHS

# flag blocks for C compilation
CFLAGS_RELEASE = -O3
CFLAGS_DEBUG   = -g -Wconversion -DDEBUG_MACROS
CFLAGS_INTERFACE  = -DPYTHON_THREADS

CFLAGS := $(CFLAGS_RELEASE)
INCLUDEPATHS = $(MULTISTRAND_INCLUDES) $(PYTHON_INCLUDES)
LIBRARYPATHS = -L/opt/local/lib


#hmm, it's not actually linking as if it were a lib, so need the full path here:

# CHANGE THIS OPTION TO REFLECT YOUR PATH TO THE libboost_python-mt shared library.
BOOSTLIB = /opt/local/lib/libboost_python-mt.a

#BOOSTLIB= /usr/lib/libboost_python-mt.so
# Above option is for Chris's system, or others where you have the .so rather than .a.
# 


LIBRARIES= $(LIBRARYPATHS) -lpython2.6
LIB_INTERFACE =  $(BOOSTLIBPATH) -shared $(BOOSTLIB)

CC = g++
COMPILE = $(CC) $(CFLAGS) $(INCLUDEPATHS)
LINK = $(CC) $(CFLAGS) $(LIBRARIES)



all: Multistrand-internal

.PHONY: all interface
# primary targets

.PHONY: debug interface-debug
# debug targets

.PHONY: clean interface-clean distclean
# cleaning targets

.PHONY: dircheck dircheck-interface Multistrand-internal 
# utilities

# .PHONE: These rules MUST run their commands.  
#
# Note that Multistrand-internal is a phony rule used to build the
# directories. If someone uses the target Multistrand, they get what
# they deserve?
# 

.SUFFIXES:
# clear out all the implicit rules that might be run.

# Targets for cleaning up our builds.
clean:
	@echo Cleaning up old object files, executables. 
	-rm -f core
	-rm -f Multistrand
	-rm -f $(OBJECTS)
	-rm -f $(MAIN_OBJECT)

interface-clean:
	@echo Cleaning up old object files, shared libraries.
	-rm -f $(INTERFACE_OBJECTS)
	-rm -f multistrand.so

distclean: interface-clean clean
	@echo Removing object file directories.
	-rmdir obj/interface/
	-rmdir obj/

# targets for setting up debugging flags.

debug: CFLAGS := $(CFLAGS_DEBUG) 
debug: clean Multistrand

interface-debug: CFLAGS := $(CFLAGS_DEBUG)
interface-debug: interface-clean interface

# utilities targets
dircheck-interface: dircheck
	@if [ ! -d obj/interface ]; then mkdir obj/interface; echo Creating obj/interface/; fi

dircheck:
	@if [ ! -d obj ]; then mkdir obj; echo Creating obj/; fi

Multistrand-internal: dircheck Multistrand

# How to build the object files.

$(INTERFACE_OBJECTS): OBJPATH := obj/python
$(INTERFACE_OBJECTS): obj/interface/%.o: %.cc $(INCLUDES)
	$(COMPILE) -c $< -o $@
# note that $< is name of 1st prereq, e.g. the .cc for which this is the .o

$(OBJPATH)/%.o: %.cc $(INCLUDES)
	$(COMPILE) $< -c -o $@


interface: INCLUDEPATHS := $(INCLUDEPATHS) $(BOOST_INCLUDES)
interface: dircheck-interface $(INTERFACE_OBJECTS)
#@echo CFLAGS=$(CFLAGS) COMPILE=$(COMPILE)
	-rm -f multistrand.so
	$(LINK) $(CFLAGS_INTERFACE) $(filter-out dircheck-interface,$^) -o multistrand.so $(LIB_INTERFACE)
#-shared $(BOOSTLIB)



Multistrand: $(MAIN_OBJECT) $(OBJECTS)
	@echo Rule expansion: $@ is target, [$?] deps newer, all: [$^].
	rm -f Multistrand
	$(LINK) $(CFLAGS) $(LIBRARIES) $(filter-out dircheck,$^) -o Multistrand

# old used: $(patsubst %,$(OBJPATH)/%,$(filter-out dircheck,$(^F))) 

# Indirectly more also depend on loop.h via headers: scomplex.h
# ssystem.h strandordering.h should generate these directly via a
# script, see note in todos.org:
# [[id:38BF8831-172D-4BC3-8B7A-D6B2EA95FE22][Code]]
#
# This is true of several others as well.

# this block currently does not work as it uses the wrong .o names.
# scomplex.o loop.o viennaenergymodel.o nupackenergymodel.o energymodel.o: energymodel.h
# loop.o move.o scomplex.o: loop.h move.h
# scomplex.o ssystem.o scomplexlist.o strandordering.o: scomplex.h optionlists.h
# scomplex.o strandordering.o: strandordering.h
# ssystem.o testingmain.o: ssystem.h
# scomplexlist.o ssystem.o: scomplexlist.h



