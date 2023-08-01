# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
Options module to grab the options object into this namespace.
"""

from ._options.options   import Options, Literals, Energy_Type
from ._options.interface import Result

Options.__module__ = 'multistrand.options'
Result.__module__ = 'multistrand.options'
