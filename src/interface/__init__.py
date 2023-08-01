# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from importlib import metadata
from subprocess import run, PIPE


__version__ = metadata.version('multistrand')

__all__ = ['objects','options','system','utils','experiment', 'concurrent', 'builder']

_call = lambda cmd, input: run(
    cmd, input=input, stdout=PIPE, stderr=PIPE, encoding="utf-8", check=True)
