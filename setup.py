# Copyright (c) 2010-2017 Caltech. All rights reserved.
# Joseph Schaeffer (schaeffer@dna.caltech.edu)
# Frits Dannenberg (fdann@dna.caltech.edu)

from distutils.core import setup, Extension
import distutils.sysconfig
import os, sys

config_vars = distutils.sysconfig.get_config_vars()

sources = ["src/system/utility.cc",
           "src/system/sequtil.cc",
           "src/interface/multistrand_module.cc",
           "src/interface/optionlists.cc",
           "src/interface/options.cc",
           "src/loop/move.cc",
           "src/loop/moveutil.cc",
           "src/loop/loop.cc",
           "src/system/energyoptions.cc",
           "src/energymodel/nupackenergymodel.cc",
           "src/energymodel/energymodel.cc",
           "src/state/scomplex.cc",
           "src/state/scomplexlist.cc",
           "src/system/statespace.cc",
           "src/system/simoptions.cc",
           "src/system/ssystem.cc",
           "src/state/strandordering.cc"
           ]

def setup_ext( ):    

    multi_ext = Extension("multistrand.system",
                          sources=sources,
                          include_dirs=["./src/include"],
                          language="c++",
                        undef_macros=['NDEBUG'],
                        extra_compile_args = ['-O3','-g', '-w', "-std=c++11", ], #FD: adding c++11 flag 
                          )
    return multi_ext



""" To update versioning, change here and in src/interface/__init__.py, and README.md """
if __name__ == '__main__':

    multi_ext = setup_ext( )
    
    setup(name="multistrand", version="2.1",
          packages=['multistrand','multistrand._options','multistrand._objects','nupack'],
          url='http://www.multistrand.org',
          license='MIT',
          author='The Multistrand Team',
          author_email='help@multistrand.org',
          package_dir={'multistrand':'src/interface'},
          ext_modules=[multi_ext])



