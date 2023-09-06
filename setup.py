# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from setuptools import setup, Extension

debug = False

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
    define_macros=[("DEBUG_MACROS", None)] if debug else [("NDEBUG", None)],
    undef_macros=["NDEBUG"] if debug else [],
    extra_compile_args=["-std=c++11", "-w", "-Wall"] + (
        ["-g", "-O0", "-fno-inline"] if debug else ["-O3"]),
    sources=[f"src/{d}/{f}.cc" for d, fs in sources.items()
             for f in fs.split(" ")])])
