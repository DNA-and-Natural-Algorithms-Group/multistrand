#
# Copyright (c) 2007-2010 Caltech. All rights reserved.
#   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)   [maintainer]
#  	             Chris Berlind (cberlind@dna.caltech.edu)
#

# A useful line for debugging dependencies:
#    @echo Rule expansion: $@ is target, [$?] deps newer, all: [$^].
#
# Note that all python build deps and so on are handled by python's distutils now, this 
#  Makefile just invokes it correctly and places the .o / built version in known locations.
# 

vpath %.h include
INCLUDES = energymodel.h loop.h move.h optionlists.h options.h python_options.h scomplex.h scomplexlist.h ssystem.h strandordering.h

SOURCES_LOOP = loop.cc move.cc
SOURCES_ENERGYMODEL = energymodel.cc nupackenergymodel.cc viennaenergymodel.cc
SOURCES_STATE = scomplex.cc scomplexlist.cc strandordering.cc
SOURCES_SYSTEM = ssystem.cc
SOURCES_OPTIONS = optionlists.cc python_options.cc
SOURCES_TESTING = testingmain.cc

SOURCES = $(SOURCES_LOOP) $(SOURCES_ENERGYMODEL) $(SOURCES_STATE) $(SOURCES_SYSTEM) $(SOURCES_OPTIONS) $(SOURCES_TESTING)

# separate the objects to make it a bit easier on ourselves when compiling targets
#  for the python library package.

VPATH=loop state system energymodel include interface

MAIN_OBJECT = testingmain.o
OBJECTS = loop.o scomplex.o energymodel.o viennaenergymodel.o nupackenergymodel.o move.o ssystem.o scomplexlist.o strandordering.o python_options.o optionlists.o


OBJPATH=obj
OBJECTS := $(OBJECTS:%=$(OBJPATH)/%)
MAIN_OBJECT := $(MAIN_OBJECT:%=$(OBJPATH)/%)

MULTISTRAND_INCLUDES := -I $(realpath ./include)

PYTHON_INCLUDES 	  = -I /usr/include/python2.6

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
PYTHON_NAMES       := python
PYTHON_DEBUG_NAMES := python2.6-dbg python-dbg
PYTHON_COMMAND = $(call find_executable,$(PYTHON_NAMES))
PYTHON_DEBUG_COMMAND = $(call find_executable,$(PYTHON_DEBUG_NAMES) $(PYTHON_NAMES))

ifeq ($(PYTHON_COMMAND),)
$(error Could not find any python executable in your PATH.)
endif

ifeq ($(call find_executable,$(PYTHON_DEBUG_NAMES)),)
$(warning Could not find a debugging python executable in your PATH. Compiling the extension module with a standard python executable.)
endif

# flag blocks for C compilation
CFLAGS_RELEASE = -O3
CFLAGS_DEBUG   = -g -Wconversion -DDEBUG_MACROS -DDEBUG
#CFLAGS_INTERFACE  = -DPYTHON_THREADS 

CFLAGS := $(CFLAGS_RELEASE)
INCLUDEPATHS = $(MULTISTRAND_INCLUDES) $(PYTHON_INCLUDES)

LIBRARIES= $(LIBRARYPATHS) -lpython2.6

# note that LIB_INTERFACE stuff is now all handled by python's distutils, we don't need to
# worry about any boost stuff, etc.
#LIB_INTERFACE = -shared $(BOOSTLIB)

CC = g++
COMPILE = $(CC) $(CFLAGS) $(INCLUDEPATHS)
LINK = $(CC) $(CFLAGS) $(LIBRARIES)



all: Multistrand-internal

.PHONY: all package
# primary targets

.PHONY: debug package-debug
# debug targets

.PHONY: clean package-clean package-debug-clean distclean
# cleaning targets

.PHONY: dircheck Multistrand-internal 
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

package-clean:
	@echo Cleaning up old object files, shared libraries.
	$(PYTHON_COMMAND) setup.py clean -b ./ -t obj/package/ --build-lib ./
# NOTE: DO NOT USE --all IN THE ABOVE COMMAND! Perhaps later if we build binaries to ./bin
# and libraries to ./lib it'll be possible, but right now that may
# delete your distribution.
	-rm -rf multistrand/

package-debug-clean:
	@echo Cleaning up old object files, shared libraries.
	$(PYTHON_DEBUG_COMMAND) setup.py clean -b ./ -t obj/package_debug/ --build-lib ./
# NOTE: DO NOT USE --all IN THE ABOVE COMMAND! Perhaps later if we build binaries to ./bin
# and libraries to ./lib it'll be possible, but right now that may
# delete your distribution.
	-rm -rf multistrand/

distclean: package-clean package-debug-clean clean
	@echo Removing object file directories.
	-rmdir obj/package_debug/
	-rmdir obj/package/
	-rmdir obj/

# Package build targets
package:
	@echo Building the 'multistrand' Python package.
	@if [ -d obj/package_debug/ ]; then $(MAKE) package-debug-clean; fi
	$(PYTHON_COMMAND) setup.py build -b ./ -t obj/package/ --build-lib ./
	@echo Package is now [hopefully] built, you can import it via "import multistrand" if the current directory is in your sys.path. In the future you may be able to run 'make install' to have it installed in your Python site packages.

package-debug:
	@echo Building the 'multistrand' Python package, with debugging symbols enabled.
	@if [ -d obj/package/ ]; then $(MAKE) package-clean; fi
	$(PYTHON_DEBUG_COMMAND) setup.py build -b ./ -t obj/package_debug/ --build-lib ./ --debug
	@echo Package is now [hopefully] built, you can import it via "import multistrand" if the current directory is in your sys.path. In the future you may be able to run 'make install' to have it installed in your Python site packages.

# targets for setting up debugging flags.

debug: CFLAGS := $(CFLAGS_DEBUG) 
debug: clean Multistrand

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

# Indirectly more also depend on loop.h via headers: scomplex.h
# ssystem.h strandordering.h should generate these directly via a
# script, see note in todos.org:
# [[id:38BF8831-172D-4BC3-8B7A-D6B2EA95FE22][Code]]
#
# This is true of several others as well.
