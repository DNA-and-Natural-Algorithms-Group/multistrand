# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from typing import Dict, Tuple

import numpy as np
import pytest

from multistrand.objects import Strand, Complex, Domain
from multistrand.options import Options
from multistrand.concurrent import MergeSim
from multistrand.utils import GAS_CONSTANT, C2K
from nupack import pfunc


@pytest.mark.parametrize("seq", ["GGGGAAACCCC"])
@pytest.mark.parametrize("celsius", [25.0, 37.0])
@pytest.mark.parametrize("rate_model", ["DNA23Metropolis", "DNA23Arrhenius"])
class Test_Equilibrium:
    @classmethod
    def test_partition_function(cls, seq: str, celsius: float, rate_model: str):
        """
        This test compares the "long-run" partition function obtained from
        Multistrand with the partition function from NUPACK.

        Note that the Multistrand estimate here is constructed by summing over
        unique states found as final states in simulated trajectories. In
        particular, this procedure discards the frequency information and does
        not constitute a Monte Carlo estimate.
        """
        kelvin = celsius + C2K
        RT = GAS_CONSTANT * kelvin

        eqd = cls.equilibrium_distribution(seq, kelvin, RT, rate_model)
        Z = 0.0
        for structure, (energy, count) in eqd.items():
            print(f"{structure:s}: {energy: e} [{count:3d} x]")
            Z += np.exp(- energy / RT)

        print()
        F = -RT * np.log(Z)
        print(f"Partition function from Multistrand trajectories: {F:e}")
        pf = pfunc([seq], material="dna", T=celsius)
        print(f"Partition function from NUPACK dynamic program  : {pf:e}")
        assert np.allclose(pf, F, rtol=1e-2)

    @classmethod
    def equilibrium_distribution(cls, seq: str, T: float, RT: float,
                                 rate_model: str) -> Dict[str, Tuple[float, int]]:
        num_sims = 1000
        sim_time = 1e-3
        sim = MergeSim()
        sim.setNumOfThreads(5)
        sim.setOptionsFactory6(cls.setup_options, num_sims, seq, T,
                               2 * np.exp(-5.06 / RT), sim_time, rate_model)
        sim.setPassageMode()
        sim.run()

        hist = {}
        for end_state in sim.results.endStates:
            for (_, _, _, _, structure, energy, _) in end_state:
                if structure in hist:
                    hist[structure][1] += 1
                else:
                    hist[structure] = [energy, 1]
        return {s: tuple(d) for s, d
                in sorted(hist.items(), key=lambda s: s[1][1], reverse=True)}

    @staticmethod
    def setup_options(num_sims: int, seq: str, temp: float, conc: float,
                      sim_time: float, rate_model: str) -> Options:
        """
        Creates an `Options` object using the sequence passed as a single domain
        with initially unpaired structure.
        """
        d = Domain(name="initial", sequence=seq, length=len(seq))
        s = Strand(domains=[d])
        c = Complex(strands=[s], structure=".")
        o = Options(
            simulation_mode="Normal", parameter_type="Nupack", substrate_type="DNA",
            num_sims=num_sims, sim_time=sim_time, start_state=[c],
            temperature=temp, join_concentration=conc)
        getattr(o, rate_model)()
        return o


if __name__ == "__main__":
    Test_Equilibrium.test_partition_function(
        "GGGGAAACCCC", 37.0, "DNA23Metropolis")
