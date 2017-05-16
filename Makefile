# Copyright (c) 2007-2017 Caltech. All rights reserved.
#   
#	Joseph Schaeffer (schaeffer@dna.caltech.edu)  
#  	Chris Berlind (cberlind@dna.caltech.edu)
#	Frits Dannenberg (fdann@caltech.edu)
#
#  All python build dependencies are handled by python distutils. 
#  This makefile invokes it and copies the built version.

vpath %.h include
INCLUDES = utility.h sequtil.h energymodel.h loop.h move.h moveutil.h optionlists.h python_options.h scomplex.h scomplexlist.h energyoptions.h simoptions.h ssystem.h strandordering.h 

SOURCES_LOOP = loop.cc move.cc moveutil.cc
SOURCES_ENERGYMODEL = energymodel.cc nupackenergymodel.cc 
SOURCES_STATE = scomplex.cc scomplexlist.cc strandordering.cc
SOURCES_SYSTEM =  utility.cc sequtil.cc  energyoptions.cc simoptions.cc ssystem.cc 
SOURCES_OPTIONS = optionlists.cc python_options.cc
SOURCES_TESTING = testingmain.cc

SOURCES = $(SOURCES_LOOP) $(SOURCES_ENERGYMODEL) $(SOURCES_STATE) $(SOURCES_SYSTEM) $(SOURCES_OPTIONS) $(SOURCES_TESTING)

VPATH=loop state system energymodel include interface
MAIN_OBJECT = testingmain.o
OBJECTS = sequtil.o loop.o scomplex.o energymodel.o  nupackenergymodel.o move.o moveutil.o energyoptions.o  simoptions.o ssystem.o scomplexlist.o strandordering.o python_options.o optionlists.o 

OBJPATH = obj
OBJECTS := $(OBJECTS:%=$(OBJPATH)/%)
MAIN_OBJECT := $(MAIN_OBJECT:%=$(OBJPATH)/%)

MULTISTRAND_INCLUDES := -I $(realpath ./include)
PYTHON_INCLUDES 	  = -I /usr/include/python2.7

## utility function for searching a path: 
## [main idea via the make info manual] 
## 1st param is file to search for, 2nd is space separated list
##  of paths (e.g. PATH environment variable with :'s changed to spaces)
## Works by clever use of the wildcard function's expansion to
## all matching files. addsuffix adds it to each path, so wildcard checks
## all the combinations and firstword is the one [in standard PATH precedence]
## that'd take priority, if it exists.

any_path_search = $(firstword $(wildcard $(addsuffix /$(1), $(2))))
path_search = $(call any_path_search,$(1),$(subst :, ,$(PATH)))

## use via, eg PYTHON := $(call path_search,python)
## or          PYTHON := $(call any_path_search,/opt/local/bin /usr/bin)

find_executable = $(strip $(firstword $(foreach file,$(1),$(call path_search,$(file)))))

## The next two blocks try to work out either a usable python version
## or if it can't find that it will just use your base python version.
##
## The debug symbols compiled version it tries some default names,
## which should be in the path. If you have a differently named version,
## or one that's not in the path, consider sym linking it in ~/bin or some
## other location in your path. E.g., I did the following:
##
## ln -s ~/Python2.6_debug/bin/python2.6 ~/bin/python2.6-dbg
## 
## Which places it in my path and with one of the searched-for names.
## If there's a better 'default' set of names for the debug version,
## we can add it to the list checked.
##
PYTHON_NAMES       := python python2.7
PYTHON_DEBUG_NAMES := python-debug
PYTHON_COMMAND = $(call find_executable,$(PYTHON_NAMES))
PYTHON_DEBUG_COMMAND = $(call find_executable,$(PYTHON_DEBUG_NAMES) $(PYTHON_NAMES))

ifeq ($(PYTHON_COMMAND),)
$(error Could not find any python executable in your PATH.)
endif

ifeq ($(call find_executable,$(PYTHON_DEBUG_NAMES)),)
$(warning Could not find a debugging python executable in your PATH. Compiling the extension module with a standard python executable.)
endif

# flag blocks for C compilation
CFLAGS_RELEASE = -O3 -std=c++11
CFLAGS_DEBUG   = -g -Wconversion -DDEBUG_MACROS -DDEBUG

#CFLAGS := $(CFLAGS_DEBUG)
CFLAGS := $(CFLAGS_RELEASE)
INCLUDEPATHS = $(MULTISTRAND_INCLUDES) $(PYTHON_INCLUDES)

LIBRARIES= $(LIBRARYPATHS) -lpython2.7

CC = g++
#CC = clang
COMPILE = $(CC) $(CFLAGS) $(INCLUDEPATHS)
LINK = $(CC) $(CFLAGS) $(LIBRARIES)


all: package

Multistrand: Multistrand-internal
# Multistrand-internal is a dummy rule to build directories. 

.PHONY: all package Multistrand
# primary targets

.PHONY: debug package-debug
# debug targets

.PHONY: clean package-clean package-debug-clean distclean
# cleaning targets

.PHONY: dircheck Multistrand-internal 
# utilities

.PHONY: docs
# documentation

.SUFFIXES:
# clear out all the implicit rules that might be run.

# Targets for cleaning up our builds.
clean: package-clean

package-clean:
	@echo Cleaning up old object files, shared libraries.
	$(PYTHON_COMMAND) setup.py clean -b ./ -t obj/package/ --build-lib ./
	-rm -rf multistrand/
	# Do not use --all here! This could delete your distribution.

distclean: package-clean package-debug-clean clean
	@echo Removing object file directories.
	-rmdir obj/package_debug/
	-rmdir obj/package/
	-rmdir obj/

# Package build targets
package:
	@echo Building the 'multistrand' Python package.
	@if [ -d obj/package_debug/ ]; then $(MAKE) package-debug-clean; fi
	@if [ -d obj/package_profiler/ ]; then $(MAKE) package-profiler-clean; fi
	$(PYTHON_COMMAND) setup.py build -b ./ -t obj/package/ --build-lib ./ --debug
	@echo Multistrand is now built. Run 'sudo make install' to install Multistrand to your Python site packages.

#documentation
docs:
	@cd doc/ && $(MAKE) clean; $(MAKE) html

#install
install:
	@echo Installing the 'multistrand' Python package to your python site-packages.
	$(PYTHON_COMMAND) setup.py install


# utilities targets
dircheck:
	@if [ ! -d obj ]; then mkdir obj; echo Creating obj/; fi

Multistrand-internal: dircheck Multistrand

# How to build the object files.
$(OBJPATH)/%.o: %.cc $(INCLUDES) 
	$(COMPILE) $< -c -o $@

Multistrand: $(MAIN_OBJECT) $(OBJECTS)
	rm -f Multistrand
	$(LINK) $(CFLAGS) $(LIBRARIES) $(filter-out dircheck,$^) -o Multistrand
