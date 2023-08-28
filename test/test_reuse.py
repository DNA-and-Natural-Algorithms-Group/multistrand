# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from contextlib import nullcontext
from typing import List, Tuple

import numpy as np
import pytest

from multistrand.objects import Complex, Domain, Strand
from multistrand.options import Options, Literals, Energy_Type
from multistrand.system import initialize_energy_model, energy
from multistrand.objects import Complex, Strand

if __name__ == "__main__":
    opt = Options(verbosity=0)
    opt.DNA23Metropolis()

    c = Complex(strands=[Strand(name="hairpin", sequence="ACT")], structure="...")
    energy = energy([c], opt, Energy_Type.Complex_energy)
    print(energy)

