#python2.6 interface/setup.py build -n -b ./ -t obj/interface/ --build-lib ./ --debug --force
# Click <mouse-2> on a history reference to select it.
# In this buffer, type RET to select the history reference near point.

# Possible history references are:
# python2.6 interface/setup.py build -n -b ./ -t obj/interface/ --build-lib ./ --debug --force
# python2.6 interface/setup.py build -n -b . -t obj/interface/ --debug --force
# python2.6 interface/setup.py build --help
# python2.6 interface/setup.py build -n -b bin/ --debug --force
# python2.6 interface/setup.py build -n -b bin/--debug --force
# python2.6 interface/setup.py build --help
# python2.6 interface/setup.py install --help
# python2.6 interface/setup.py install -O0 -n --root=./ --install-lib=. --prefix=None
# python2.6 interface/setup.py install -O0 -n --root=./ --install-lib=. --prefix
# python2.6 interface/setup.py install -O0 -n --root=./ --install-lib=.
# python2.6 interface/setup.py install -O0 -n --install-base=lib/ --root=./ --install-lib=.
# python2.6 interface/setup.py install -O0 -n --install-base=lib/ --install-lib=.
# python2.6 interface/setup.py install -O0 -n --install-base=lib/
# ls bin/lib/python2.6/site-packages/
# python2.6 interface/setup.py install -O0 -n --root=./ --prefix=bin/
# python2.6 interface/setup.py install -O0 -n --root=./
# python2.6 interface/setup.py install --help
# python2.6 setup.py install --help
# ln -s build/lib.macosx-10.6-x86_64-2.6/multistrand.so ./multistrand.so
# rm multistrand.so
# python2.6 interface/setup.py build --debug --force
# man gcc
# python2.6 interface/setup.py build --debug --force
# echo $CFLAGS 	echo $DISTUTILS_DEBUG
# python2.6 interface/setup.py build --debug --force
# export DISTUTILS_DEBUG=1
# make interface-debug 	make interfare-debug
# grep Makefile
# python2.6 interface/setup.py build --debug --force
# python2.6 interface/setup.py build --help-compiler
# python2.6 interface/setup.py build --help
# python2.6 interface/setup.py build
# python2.6 interface/setup.py clean --all
# python2.6 interface/setup.py clean --help
# python2.6 interface/setup.py --help-commands
# python2.6 interface/setup.py cmd --help
# python2.6 interface/setup.py distclean
# python2.6 interface/setup.py build
# python2.6 interface/setup.py clean
# python2.6 interface/setup.py build
# ls -la /usr/bin/g++* 	ls /usr/bin/g++*
# which g++
# ln -s build/lib.macosx-10.6-x86_64-2.6/multistrand.so ./multistrand.so
# python2.6 interface/setup.py build
# setup.py for multistrand module.
from distutils.core import setup, Extension

multi_ext = Extension("multistrand",
                      sources=["interface/multistrand_module.cc",
                               "loop/move.cc",
                               "loop/loop.cc",
                               "energymodel/nupackenergymodel.cc",
                               "energymodel/energymodel.cc",
                               "state/scomplex.cc",
                               "state/scomplexlist.cc",
                               "system/ssystem.cc",
                               "state/strandordering.cc",
                               "energymodel/viennaenergymodel.cc",
                               "optionlists.cc",
                               "python_options.cc"
                               ],
                      include_dirs=["./include"],
                      language="c++",
                      define_macros=[('DEBUG',None),
                                     ('DEBUG_MACROS',None)],
                      undef_macros=['NDEBUG'],
                      extra_compile_args = ['-w'],
                      )

setup(name="multistrand", version="1.0",
      ext_modules=[multi_ext])
