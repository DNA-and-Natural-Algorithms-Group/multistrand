#
# Copyright (c) 2007-2010 Caltech. All rights reserved.
#   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)   [maintainer]
#  	             Chris Berlind (cberlind@dna.caltech.edu)
#

vpath %.h include
INCLUDES = energymodel.h loop.h move.h optionlists.h options.h python_options.h scomplex.h scomplexlist.h ssystem.h strandordering.h

SOURCES_LOOP = loop/loop.cc loop/move.cc
SOURCES_ENERGYMODEL = energymodel/energymodel.cc energymodel/nupackenergymodel.cc energymodel/viennaenergymodel.cc
SOURCES_STATE = state/scomplex.cc state/scomplexlist.cc state/strandordering.cc
SOURCES_SYSTEM = system/ssystem.cc
SOURCES_OPTIONS = optionlists.cc python_options.cc
SOURCES_TESTING = python_testing/embedding_test.cc python_testing/embedding_test2.cc testingmain.cc

SOURCES = $(SOURCES_LOOP) $(SOURCES_ENERGYMODEL) $(SOURCES_STATE) $(SOURCES_SYSTEM) $(SOURCES_OPTIONS) $(SOURCES_TESTING)

# separate the objects to make it a bit easier on ourselves when compiling targets
#  for the python library interface.

VPATH=obj loop state system energymodel

MAIN_OBJECT = testingmain.o
OBJECTS = loop.o scomplex.o energymodel.o viennaenergymodel.o nupackenergymodel.o move.o ssystem.o scomplexlist.o strandordering.o python_options.o optionlists.o
PYTHON_OBJECTS := $(OBJECTS:%=obj/python/%)
OBJPATH=obj

INCLUDEFLAGS = -I ./include
CFLAGS = -O3 -Wconversion $(INCLUDEFLAGS)
CFLAGS_DEBUG = -g -Wconversion $(INCLUDEFLAGS)

# flags used for the python interface compilation only
CFLAGS_PYTHON = -O3 -DPYTHON_THREADS $(INCLUDEFLAGS) -I /usr/include/python2.6

LIBRARIES= -lpython2.6

CC = g++
COMPILE = $(CC) $(CFLAGS)



all: Multistrand

.PHONY: all clean python-interface debug debug-full distclean FORCE

FORCE:

python-interface: FORCE $(PYTHON_OBJECTS) 

distclean: clean python-clean

clean: 
	-rm *.o core
	-rm iosys/*.o

python-clean:
	-rm $(PYTHON_OBJECTS)
	-rmdir python-interface

debug: FORCE
	$(MAKE) clean
	$(MAKE) debug-full


debug-full: CFLAGS := $(CFLAGS_DEBUG) 
debug-full: Multistrand

$(OBJECTS): obj/.
obj/.: 
	mkdir obj

obj/python/.: obj/.
	mkdir obj/python

$(PYTHON_OBJECTS): OBJPATH := obj/python
$(PYTHON_OBJECTS): CFLAGS :=  $(CFLAGS_PYTHON)
$(PYTHON_OBJECTS): obj/python/%: %
# note that if this rule is checked, any dependency newer than the python .o file will be recalculaed.

#%.o: FORCE
#	$(CC) $(CFLAGS_PYTHON) -c $(patsubst obj/python/%.o,%.cc,$@) -o $@

Multistrand: $(MAIN_OBJECT) $(OBJECTS)
	rm -f Multistrand
	$(CC) $(CFLAGS) $(LIBRARIES) $(patsubst %,$(OBJPATH)/%,$(^F)) -o Multistrand

# Indirectly more also depend on loop.h via headers: scomplex.h
# ssystem.h strandordering.h should generate these directly via a
# script, see note in todos.org:
# [[id:38BF8831-172D-4BC3-8B7A-D6B2EA95FE22][Code]]
#
# This is true of several others as well.

scomplex.o loop.o viennaenergymodel.o nupackenergymodel.o energymodel.o: energymodel.h
loop.o move.o scomplex.o: loop.h move.h
scomplex.o ssystem.o scomplexlist.o strandordering.o: scomplex.h optionlists.h
scomplex.o strandordering.o: strandordering.h
ssystem.o testingmain.o: ssystem.h
scomplexlist.o ssystem.o: scomplexlist.h

%.o: obj/.
%.o: %.cc
	$(COMPILE) -c $< -o $(OBJPATH)/$@

