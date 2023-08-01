# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from itertools import chain
from typing import List

import numpy as np
import pytest

from multistrand.objects import Complex, Domain, Strand
from multistrand.options import Options, Literals
from multistrand.system import SimSystem


class Test_Determinism:
    """
    Compare two Multistrand trajectories which are configured deterministically.
    """
    # jump indices to compare between trajectories
    history = np.s_[-int(1e4)::int(1e3)]

    @pytest.mark.parametrize("seed", list(np.random.randint(int(1e10), size=10)))
    @pytest.mark.parametrize(
        "toehold_seq, bm_design_B, toehold_extra, structure",
        [("GTGGGT", "ACCGCACGTCACTCACCTCG", "TTT",
          "..(.((.....)).).....((((((+))))))((((((((((((((((((((+))))))))))))))))))))")
         ])
    def test_determinism(cls, toehold_seq: str, bm_design_B: str, toehold_extra: str,
                         structure: str, seed: int) -> None:
        # build complexes with domain-level information
        toehold = Domain(
            name="toehold", sequence=toehold_seq, length=len(toehold_seq))
        branch_migration_B = Domain(
            name="bm_B", sequence=bm_design_B, seq_length=len(bm_design_B))
        substrate_B = toehold + branch_migration_B
        incumbent_B = Strand(name="incumbent", domains=[branch_migration_B.C])
        incoming_B = substrate_B.C
        start_complex_B1 = Complex(
            strands=[incoming_B, substrate_B, incumbent_B], structure=structure)

        # compare trajectories
        assert (cls.single_run(start_complex_B1, seed, True)
                == cls.single_run(start_complex_B1, seed, False))

    @classmethod
    def single_run(cls, start_complex: Complex,
                   seed: int, verbose: bool) -> List[float | str]:
        opt = cls.create_config(start_complex, seed)
        sys = SimSystem(opt)
        sys.start()
        return cls.summarise_trajectory(opt, verbose)

    @staticmethod
    def create_config(start_complex: Complex, seed: int) -> Options:

        o2 = Options()
        o2.simulation_mode = Literals.trajectory
        o2.initial_seed = seed
        o2.num_simulations = 1
        o2.simulation_time = 0.0001
        o2.temperature = 37.0
        o2.dangles = 1
        o2.start_state = [start_complex]
        o2.output_interval = 1

        o2.JSDefault()
        return o2

    @classmethod
    def summarise_trajectory(cls, opt: Options, verbose: bool) -> List[float | str]:
        seqstring = ''
        times = opt.full_trajectory_times
        tubestructs = []
        energies = []

        # go through each output microstate of the trajectory
        for i in range(len(opt.full_trajectory)):
            # this is a list of the complexes present in this tube microstate
            states = opt.full_trajectory[i]
            # similarly, extract the secondary structures for each complex
            structs = [s[4] for s in states]
            # give the dot-paren secondary structure for the whole test tube
            tubestructs.append(' '.join(structs))

            dG = sum((s[5] for s in states))
            energies.append(dG)

            # extract the strand sequences in each complex
            # (joined by "+" for multistranded complexes)
            newseqs = [s[3] for s in states]
            # make a space-separated string of complexes, to represent the whole tube system sequence
            newseqstring = ' '.join(newseqs)

            if not newseqstring == seqstring:
                # because strand order can change upon association or
                # dissociation, print it when it changes
                seqstring = newseqstring

        summary = list(chain(times[cls.history], tubestructs[cls.history],
                             energies[cls.history]))
        if verbose:
            [print(x) for x in summary]
            print()
        return summary


# ==============================================================================


if __name__ == "__main__":
    toehold_seq = "GTGGGT"
    bm_design_B = "ACCGCACGTCACTCACCTCG"
    toehold_extra = "TTT"
    structure = "..(.((.....)).).....((((((+))))))((((((((((((((((((((+))))))))))))))))))))"
    seed = np.random.randint(int(1e10))
    Test_Determinism().test_determinism(
        toehold_seq, bm_design_B, toehold_extra, structure, seed)
