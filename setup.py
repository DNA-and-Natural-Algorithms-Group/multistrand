# Copyright (c) 2010-2017 Caltech. All rights reserved.
# Joseph Schaeffer (schaeffer@dna.caltech.edu)
# Frits Dannenberg (fdann@dna.caltech.edu)

from distutils.core import setup, Extension
import distutils.sysconfig
import os, sys

config_vars = distutils.sysconfig.get_config_vars()

sources = ["system/utility.cc",
           "system/sequtil.cc",
           "interface/multistrand_module.cc",
           "interface/optionlists.cc",
           "interface/options.cc",
           "loop/move.cc",
           "loop/moveutil.cc",
           "loop/loop.cc",
           "system/energyoptions.cc",
           "energymodel/nupackenergymodel.cc",
           "energymodel/energymodel.cc",
           "state/scomplex.cc",
           "state/scomplexlist.cc",
           "system/simoptions.cc",
           "system/ssystem.cc",
           "state/strandordering.cc"
           ]

def setup_ext( ):    

    multi_ext = Extension("multistrand.system",
                          sources=sources,
                          include_dirs=["./include"],
                          language="c++",
                        undef_macros=['NDEBUG'],
                        extra_compile_args = ['-O3','-g', '-w', "-std=c++11", ], #FD: adding c11 flag
                          )
    return multi_ext

if __name__ == '__main__':

    multi_ext = setup_ext( )
    
    setup(name="multistrand", version="2.1",
          packages=['multistrand','multistrand._options','multistrand._objects','nupack'],
          url='http://www.multistrand.org',
          license='MIT',
          author='The Multistrand Team',
          author_email='help@multistrand.org',
          package_dir={'multistrand':'interface'},
          ext_modules=[multi_ext])



