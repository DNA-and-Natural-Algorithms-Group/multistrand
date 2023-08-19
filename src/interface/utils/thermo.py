# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

import numpy as np

from .._objects.strand import Strand
import math
from nupack.analysis import pfunc as _nu_pfunc
from nupack.analysis import pairs, mfe, structure_probability, ensemble_size, subopt, sample
from nupack.analysis import energy as structure_energy
from nupack import Model as _NU_Model

# physical constants

GAS_CONSTANT = 0.0019872036
""" kcal / (K * mol) """

C2K = 273.15
""" Celsius to Kelvin """


NUPACK3 = "-nupack3"
RNA_NUPACK = "rna06" + NUPACK3
DNA_NUPACK = "dna04" + NUPACK3

class Model(_NU_Model):
    def __init__(self, option=None, ensemble="some", material="DNA", kelvin=None, celsius=37, sodium=1.0, magnesium=0.0):
        kwargs = {}
        if option is None:
            kwargs["ensemble"] = ensemble
            kwargs["material"] = material
            kwargs["kelvin"] = kelvin
            kwargs["celsius"] = celsius
            kwargs["sodium"] = sodium
            kwargs["magnesium"] = magnesium
        else:
            kwargs["ensemble"] = option.dangleToString[option.dangles]
            kwargs["material"] = option.substrateToString[option.substrate_type]
            kwargs["kelvin"] = option.temperature
            kwargs["sodium"] = option.sodium
            kwargs["magnesium"] = option.magnesium

        kwargs["material"] = RNA_NUPACK if kwargs["material"] == "RNA" else DNA_NUPACK
        kwargs["ensemble"] = kwargs["ensemble"].lower() + NUPACK3

        super(Model, self).__init__(**kwargs)

    @_NU_Model.ensemble.setter
    def ensemble(self, value):
        self._ensemble = value.lower() + NUPACK3

    @_NU_Model.material.setter
    def material(self, value):
        if value == "RNA":
            self._material = RNA_NUPACK
        else:
            self._material = DNA_NUPACK


def _dGadjust(K, N):
    """Adjust NUPACK's native free energy (with reference to mole fraction units) to be appropriate for molar units, assuming N strands in the complex."""
    water = 55.14  # molar concentration of water at 37 C, ignore temperature dependence, which is about 5%
    adjust = GAS_CONSTANT * K * math.log(water)  # converts from NUPACK mole fraction units to molar units, per association
    return adjust * (N - 1)


def pfunc(strands, model=None):
    """Overriding NUPACK's pfunc so that it returns an adjusted dG value"""
    if model is None:
        model = Model()
    return _nu_pfunc(strands=strands, model=model)[1] + _dGadjust(model.temperature, len(strands))

def meltingTemperature(seq, concentration=1.0e-9):
    """
    Returns the melting temperature in Kelvin for a duplex of the given sequence.
    Sequences should be at least 8 nt long for the SantaLucia model to reasonably apply.
    For shorter sequences, see notes on "Melting Temperature (Tm) Calculation" by biophp.org
    Specifically, look at basicTm vs Base-stacking Tm. FD mar 2018
    """
    strand = Strand(sequence=seq)
    model20 = Model(material='dna', celsius=20)
    energy20 = (float(mfe([strand.sequence, strand.C.sequence], model20)[0][1])
                + GAS_CONSTANT * (20 + C2K) * np.log(55.5))
    model30 = Model(material='dna', celsius=20)
    energy30 = (float(mfe([strand.sequence, strand.C.sequence], model30)[0][1])
                + GAS_CONSTANT * (30 + C2K) * np.log(55.5))

    dS = (energy20 - energy30) / 10.0  # kcal/ K mol
    dH = energy30 + (30.0 + C2K) * dS  # kcal/mol
    return (dH / (dS + GAS_CONSTANT * np.log(concentration / 4.0)))


def dGC_feature(o, i: int):
    states = o.full_trajectory[i]
    dGC = 0.0
    for state in states:
        dGC += (state[5] - (
                o._temperature_kelvin * 0.0019872036
                * np.log(1.0 / o.join_concentration) * state[4].count("+")))
        # print("count is  " +  str(state[4].count("+")))
        # print("join conc is " + str(o.join_concentration))
        # print("dG-Complex is " + "%.2f" % dGC + " kcal/mol  for " + str(state[3]))
    return f"dGC={dGC:> 6.2f} kcal/mol"
