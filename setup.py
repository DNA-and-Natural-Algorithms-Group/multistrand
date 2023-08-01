# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from setuptools import setup, Extension

sources = {
    "system": "ssystem simoptions energyoptions statespace utility sequtil simtimer",
    "state": "scomplex scomplexlist strandordering",
    "loop": "move moveutil loop",
    "energymodel": "nupackenergymodel energymodel",
    "interface": "multistrand_module options optionlists"
}

setup(ext_modules=[Extension(
    name="multistrand.system",
    include_dirs=["./src/include"],
    language="c++",
    undef_macros=['NDEBUG'],
    # define_macros=[('DEBUG_MACROS', None)],
    extra_compile_args = ["-O3", "-std=c++11", "-Wall", "-Wconversion"],
    sources=[f"src/{d}/{f}.cc" for d, fs in sources.items()
             for f in fs.split(" ")])])
