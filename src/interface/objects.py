# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
Objects module to get the useful names into the correct namespace.
"""

__module__ = 'multistrand.objects'

from ._objects.domain  import Domain
from ._objects.complex import Complex
from ._objects.strand import Strand
from ._objects.stopcondition import StopCondition, Macrostate

__all__ = ['Strand','StopCondition','Macrostate','Complex','Domain']

# The following appears to be necessary [at the moment] as otherwise
# we can't generate appropriate documentation for these sub objects as
# they look like they are in a different package. Perhaps we should
# change the directory hierarchy - rename _objects to objects, and
# then moving this file into the new directory as __init__.py and
# importing from local space.

Strand.__module__ = 'multistrand.objects'
Complex.__module__ = 'multistrand.objects'
Macrostate.__module__ = 'multistrand.objects'
StopCondition.__module__ = 'multistrand.objects'
Domain.__module__ = 'multistrand.objects'
