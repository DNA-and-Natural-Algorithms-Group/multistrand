#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
#
#
# setup the python distribution
#
#

sources = ["system/utility.cc",
           "interface/multistrand_module.cc",
           "interface/optionlists.cc",
           "interface/options.cc",
           "loop/move.cc",
           "loop/loop.cc",
           "energymodel/nupackenergymodel.cc",
           "energymodel/energymodel.cc",
           "state/scomplex.cc",
           "state/scomplexlist.cc",
           "system/energyoptions.cc",
           "system/simoptions.cc",
           "system/ssystem.cc",
           "state/strandordering.cc",
           "energymodel/viennaenergymodel.cc"
           ]

# setup.py for multistrand module.
from distutils.core import setup, Extension

import distutils.sysconfig
config_vars = distutils.sysconfig.get_config_vars()
# config variables used in the next two functions.


def setup_flags():
    # remove -Wstrict-prototypes from our default flags for the compiler, if it's in there
    # 
    opt = config_vars['OPT'].split()
    if '-Wstrict-prototypes' in opt:
        opt.remove('-Wstrict-prototypes')
        config_vars['OPT'] = ' '.join( opt )
    cflags = config_vars['CFLAGS'].split()
    if '-Wstrict-prototypes' in cflags:
        cflags.remove('-Wstrict-prototypes')
        config_vars['CFLAGS'] = ' '.join( cflags )
    cppflags = config_vars['CPPFLAGS'].split()
    if '-Wstrict-prototypes' in cppflags:
        cppflags.remove('-Wstrict-prototypes')
        config_vars['CPPFLAGS'] = ' '.join( cppflags )

def setup_libcheck():
    # Find the library path [so we can check for tcmalloc].
    import os, os.path
    ldflags = config_vars['LDFLAGS'].split() + ['/usr/lib']
    for libpath in ldflags + os.environ.get('LIBRARY_PATH', '').split():
        if libpath.startswith('-L'):
            pathname = libpath[2:] # strip the -L from each entry
        else:
	     pathname = libpath
        try:
            items = iter( os.listdir(pathname) )
        except OSError:
            continue
        for item in items:
            if item.startswith('libtcmalloc'):
                return True
    return False
    

def setup_ext( have_tcmalloc):    
    import sys
    if have_tcmalloc:
        tcmalloc_optional = ['tcmalloc_minimal']
        print "Found tcmalloc, will attempt to use it."
    else:
        print "Cannot find tcmalloc, proceeding with standard memory management."
        tcmalloc_optional = []
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
                              libraries=[] + tcmalloc_optional,
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
                              libraries=['tcmalloc_and_profiler'],  #google perftools profiler.
                              extra_compile_args = ['-w','-g','-O0'],
                              )
        
    else:
        multi_ext = Extension("multistrand.system",
                              sources=sources,
                              include_dirs=["./include"],
                              language="c++",
                              libraries=[] + tcmalloc_optional,
                              #                          define_macros=[('DEBUG',None),
                              #                                         ('DEBUG_MACROS',None)],
                              #                                         ('Py_TRACE_REFS',None)],
                              undef_macros=['NDEBUG'],
                              extra_compile_args = ['-O3','-g', '-w'],
                              #                          ['-Wno-strict-prototypes','-w','-O0','-v','-fcommon', '-fno-wrapv'],   #This is 'disable all warnings compiler flag'
                              )
    return multi_ext

if __name__ == '__main__':
    setup_flags()
    tcmalloc_flag = setup_libcheck()
    multi_ext = setup_ext( have_tcmalloc=tcmalloc_flag)
    
    setup(name="multistrand", version="2.0",
          packages=['multistrand','multistrand._options','multistrand._objects','nupack'],
          url='http://www.multistrand.org',
          license='MIT',
          author='The Multistrand Team',
          author_email='help@multistrand.org',
          package_dir={'multistrand':'interface'},
          ext_modules=[multi_ext])

#          zip_safe=False)


