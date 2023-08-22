# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

import math
from functools import wraps
from typing import Callable

import numpy as np
import nupack
import nupack.analysis

from .._objects.strand import Strand


# physical constants ===========================================================


GAS_CONSTANT = 0.0019872036
""" kcal / (K * mol) """

C2K = 273.15
""" Celsius to Kelvin """


# NUPACK interface =============================================================


class Model(nupack.Model):
    """
    Adapter for the `nupack.Model` in NUPACK 4, migrating the legacy NUPACK 3
    calling conventions used in Multistrand.
    """
    NUPACK3 = "-nupack3"
    PARAM = {"rna": "rna99", "dna": "dna04"}

    def __init__(self, option=None, complex=None,
                 material="DNA", ensemble="some",
                 kelvin=None, celsius=None, sodium=1.0, magnesium=0.0):
        kwargs = {}
        if option is None and complex is None:
            kwargs["material"] = material
            kwargs["ensemble"] = ensemble
            if celsius is not None and kelvin is None:
                kwargs["celsius"] = celsius
            elif kelvin is not None and celsius is None:
                kwargs["kelvin"] = kelvin
            elif celsius is None and kelvin is None:
                kwargs["celsius"] = 37.0
            else:
                raise ValueError("Please specify either `celsius` or `kelvin`.")
            kwargs["sodium"] = sodium
            kwargs["magnesium"] = magnesium
        elif option.__class__.__name__ == "Options" and complex is None:
            kwargs["material"] = option.substrateToString[option.substrate_type]
            kwargs["ensemble"] = option.dangleToString[option.dangles]
            kwargs["kelvin"] = option.temperature
            kwargs["sodium"] = option.sodium
            kwargs["magnesium"] = option.magnesium
        elif complex.__class__.__name__ == "Complex" and option is None:
            kwargs["material"] = complex._substrate_type
            kwargs["ensemble"] = complex._dangles
            kwargs["celsius"] = complex._temperature
            kwargs["sodium"] = complex._sodium
            kwargs["magnesium"] = complex._magnesium
        else:
            raise ValueError(
                "Unrecognised constructor arguments for `utils.thermo.Model`")

        kwargs["material"] = self._migrate_material(kwargs["material"])
        kwargs["ensemble"] = self._migrate_ensemble(kwargs["ensemble"])

        super().__init__(**kwargs)

    @classmethod
    def _migrate_material(cls, material: str) -> str:
        return cls.PARAM[material.lower()] + cls.NUPACK3

    @classmethod
    def _migrate_ensemble(cls, ensemble: str) -> str:
        return ensemble.lower() + cls.NUPACK3

    @classmethod
    def migrate_model(cls, **kwargs) -> "Model":
        model = kwargs.pop("model", None)
        if model is None:
            model = cls(**kwargs)
        assert isinstance(model, cls)
        return model

    @classmethod
    def migrate_nupack_call(cls, func_name: str, n_args: int) -> Callable:
        try:
            nu_func = getattr(nupack, func_name)
        except AttributeError:
            nu_func = getattr(nupack.analysis, func_name)
        @wraps(nu_func)
        def func(*args, **kwargs):
            assert len(args) == n_args, "Unexpected calling convention"
            return nu_func(*args, model=cls.migrate_model(**kwargs))
        return func

    @classmethod
    def migrate_nupack(cls) -> None:
        nupack_funcs = [
            ("energy", 2), ("ensemble_size", 1), ("mfe", 1), ("pairs", 1),
            ("pfunc", 1), ("sample", 2), ("prob", 2),
            ("structure_probability", 2), ("subopt", 2), ("defect", 2)]
        for name, n_args in nupack_funcs:
            globals()[name] = cls.migrate_nupack_call(name, n_args)


Model.migrate_nupack()


def complex_free_energy(strands, **kwargs):
    """
    Call NUPACK's pfunc, discard the partition function value and adjust units
    for the free energy value.
    """
    model = Model.migrate_model(**kwargs)
    return pfunc(strands, model=model)[1] + _dGadjust(model.temperature, len(strands))


# Multistrand utilities ========================================================


def _dGadjust(K, N):
    """
    Adjust NUPACK's native free energy (with reference to mole fraction units)
    to be appropriate for molar units, assuming N strands in the complex.
    """
    # molar concentration of water at 37 C, ignoring temperature dependence,
    # which is about 5%
    water = 55.14
    # converts from NUPACK mole fraction units to molar units, per association
    adjust = GAS_CONSTANT * K * math.log(water)
    return adjust * (N - 1)


def meltingTemperature(seq, concentration=1.0e-9):
    """
    Returns the melting temperature in Kelvin for a duplex of the given
    sequence. Sequences should be at least 8 nt long for the SantaLucia model to
    reasonably apply. For shorter sequences, see notes on "Melting Temperature
    (Tm) Calculation" by biophp.org Specifically, look at basicTm vs
    Base-stacking Tm. FD mar 2018
    """
    strand = Strand(sequence=seq)
    energy_at = lambda T: (mfe([strand.sequence, strand.C.sequence],
                               material='dna', celsius=T)[0].energy
                           + GAS_CONSTANT * (T + C2K) * np.log(55.5))
    e20, e30 = energy_at(20), energy_at(30)
    dS = (e20 - e30) / 10.0  # kcal/ K mol
    dH = e30 + (30.0 + C2K) * dS  # kcal/mol
    return dH / (dS + GAS_CONSTANT * np.log(concentration / 4.0))


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
