#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
#
#
# Objects module to get the useful names into the correct namespace.
#

__module__ = 'multistrand.objects'

from _objects.objects import Strand, RestingState, StopCondition
from _objects.complex import Complex
from _objects.domain  import Domain

Strand.__module__ = 'multistrand.objects'
Complex.__module__ = 'multistrand.objects'
RestingState.__module__ = 'multistrand.objects'
StopCondition.__module__ = 'multistrand.objects'
Domain.__module__ = 'multistrand.objects'
