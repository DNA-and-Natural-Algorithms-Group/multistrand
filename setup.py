#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
#
#
# setup the python distribution
#
#

sources = ["interface/multistrand_module.cc",
           "interface/optionlists.cc",
           "interface/options.cc",
           "loop/move.cc",
           "loop/loop.cc",
           "energymodel/nupackenergymodel.cc",
           "energymodel/energymodel.cc",
           "state/scomplex.cc",
           "state/scomplexlist.cc",
           "system/ssystem.cc",
           "state/strandordering.cc",
           "energymodel/viennaenergymodel.cc"
           ]

# setup.py for multistrand module.
from distutils.core import setup, Extension

import sys
if '--use-debug-defs' in sys.argv:
    sys.argv.remove('--use-debug-defs')
    multi_ext = Extension("multistrand.system",
                          sources=sources,
                          include_dirs=["./include"],
                          language="c++",
                          define_macros=[('DEBUG',None),
                                         ('DEBUG_MACROS',None),
                                        ('Py_TRACE_REFS',None)],
                          undef_macros=['NDEBUG'],
                          #This is 'disable all warnings compiler flag' [possibly shouldn't be used for debug version]:
                          extra_compile_args = ['-Wno-strict-prototypes','-w','-g','-O0'],
                      )
elif '--use-profiler-defs' in sys.argv:
    sys.argv.remove('--use-profiler-defs')
    multi_ext = Extension("multistrand.system",
                          sources=sources,
                          include_dirs=["./include"],
                          language="c++",
                          define_macros=[('DEBUG',None),
                                         ('DEBUG_MACROS',None),
                                        ('PROFILING',None)],
                          undef_macros=['NDEBUG'],
                          libraries=['profiler'],  #google perftools profiler.
                          extra_compile_args = ['-w','-g','-O0'],
                      )
    
else:
    multi_ext = Extension("multistrand.system",
                          sources=sources,
                          include_dirs=["./include"],
                          language="c++",
                          libraries=['tcmalloc_minimal'],
#                          define_macros=[('DEBUG',None),
#                                         ('DEBUG_MACROS',None)],
#                                         ('Py_TRACE_REFS',None)],
                          undef_macros=['NDEBUG'],
                          extra_compile_args = ['-O3','-g', '-w'],
#                          ['-Wno-strict-prototypes','-w','-O0','-v','-fcommon', '-fno-wrapv'],   #This is 'disable all warnings compiler flag'
                      )

setup(name="multistrand", version="1.0",
      packages=['multistrand','multistrand._options','multistrand._objects'],
      package_dir={'multistrand':'interface'},
      ext_modules=[multi_ext])

#      packages=['interface','interface._options','interface._objects'],
#      package_dir={'interface':''},

