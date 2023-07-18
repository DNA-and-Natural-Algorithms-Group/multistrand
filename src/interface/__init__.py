####################################################################
#                                                                  #
#                                                                  #
#   Copyright (c) 2018 Caltech. All rights reserved.               #
#   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)         #
#             Frits Dannenberg (fdann@caltech.edu)                 #
#                                                                  #
#                                                                  #
####################################################################


from importlib import metadata
from subprocess import run, PIPE


__version__ = metadata.version('multistrand')

__all__ = ['objects','options','system','utils','experiment', 'concurrent', 'builder']

_call = lambda cmd, input: run(
    cmd, input=input, stdout=PIPE, stderr=PIPE, encoding="utf-8", check=True)
