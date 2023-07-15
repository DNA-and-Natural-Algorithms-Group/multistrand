####################################################################
#                                                                  #
#                                                                  #
#   Copyright (c) 2018 Caltech. All rights reserved.               #
#   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)         #
#             Frits Dannenberg (fdann@caltech.edu)                 #
#                                                                  #
#                                                                  #
####################################################################


from subprocess import run, PIPE

_call = lambda cmd, input: run(
    cmd, input=input, stdout=PIPE, stderr=PIPE, encoding="utf-8", check=True)


# defines what 'from multistrand import *' means.
__all__ = ['objects','options','system','utils','experiment', 'concurrent', 'builder']
__version__ = "2.1"
