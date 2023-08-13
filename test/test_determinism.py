# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from contextlib import nullcontext
from typing import List, Tuple

import numpy as np
import pytest

from multistrand.objects import Complex, Domain, Strand
from multistrand.options import Options, Literals
from multistrand.utils import generate_sequence
from multistrand.system import SimSystem
from multistrand.concurrent import MergeSim


@pytest.mark.parametrize(
    "rate_model",
    ["JSMetropolis25", "JSMetropolis37", "JSKawasaki25", "JSKawasaki37",
     "DNA23Metropolis", "DNA23Arrhenius"])
class Test_Determinism:
    # jump indices to compare between trajectories
    history = np.s_[::int(1e3)]

    @pytest.mark.parametrize(
        "toehold_seq, bm_design_B, toehold_extra, structure",
        [("GTGGGT", "ACCGCACGTCACTCACCTCG", "TTT",
          "..(.((.....)).).....((((((+))))))((((((((((((((((((((+))))))))))))))))))))")
         ])
    @pytest.mark.parametrize("seed", list(np.random.randint(int(1e10), size=10)))
    def test_trajectories(cls, rate_model: str,
                          toehold_seq: str, bm_design_B: str, toehold_extra: str,
                          structure: str, seed: int) -> None:
        """
        Compare two Multistrand trajectories which are configured identically.
        """
        # build initial state
        toehold = Domain(
            name="toehold", sequence=toehold_seq, length=len(toehold_seq))
        branch_migration_B = Domain(
            name="bm_B", sequence=bm_design_B, seq_length=len(bm_design_B))
        substrate_B = toehold + branch_migration_B
        incumbent_B = Strand(name="incumbent", domains=[branch_migration_B.C])
        incoming_B = substrate_B.C
        start_complex = Complex(
            strands=[incoming_B, substrate_B, incumbent_B], structure=structure)

        # simulate trajectories
        structs1, energies1, times1 = cls.single_run(
            rate_model, [start_complex], seed, True)
        structs2, energies2, times2 = cls.single_run(
            rate_model, [start_complex], seed, False)

        # compare trajectories
        if rate_model == "DNA23Arrhenius":
            common = np.s_[:min(len(structs1), len(structs2))]
            structs1, structs2 = structs1[common], structs2[common]
            energies1, energies2 = energies1[common], energies2[common]
            times1, times2 = times1[common], times2[common]
        assert (structs1 == structs2).all()
        assert (energies1 == energies2).all()
        if rate_model == "DNA23Arrhenius":
            assert np.allclose(times1, times2, rtol=1e-2)
        else:
            assert (times1 == times2).all()

    @classmethod
    def single_run(cls, rate_model: str, start_state: List[Complex],
                   seed: int, verbose: bool) -> Tuple[np.ndarray,...]:
        opt = cls.create_config(rate_model, start_state, seed)
        sys = SimSystem(opt)
        sys.start()
        return cls.summarise_trajectory(opt, verbose)

    @staticmethod
    def create_config(rate_model: str, start_state: List[Complex],
                      seed: int) -> Options:
        o = Options()
        o.simulation_mode = Literals.trajectory
        o.initial_seed = seed
        o.num_simulations = 1
        o.simulation_time = 5e-4
        o.temperature = 37.0
        o.dangles = 1
        o.start_state = start_state
        o.output_interval = 1
        getattr(o, rate_model)()
        return o

    @classmethod
    def summarise_trajectory(cls, opt: Options,
                             verbose: bool) -> Tuple[np.ndarray,...]:
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

        summary = (np.array(tubestructs[cls.history]),
                   np.array(energies[cls.history]),
                   np.array(times[cls.history]))
        if verbose:
            for feature, fmt in zip(summary, ["s", "e", "e"]):
                [print(f"{x:{fmt}}") for x in feature]
                print()
        return summary

    @pytest.mark.parametrize(
        "random_seq, raising_random",
        [(True, pytest.raises(ValueError)), (False, nullcontext())])
    @pytest.mark.parametrize(
        "wrong_mode, raising_mode",
        [(True, pytest.raises(TypeError)), (False, nullcontext())])
    @pytest.mark.parametrize("seqlen", [7, 23])
    @pytest.mark.parametrize("seed", list(np.random.randint(int(1e10), size=2)))
    def test_options_factory(
            self, rate_model: str, seqlen: int, seed: int,
            random_seq: bool, raising_random, wrong_mode: bool, raising_mode):
        """
        Exercise precondition checks for the user-provided `OptionsFactory()`,
        whose output should be deterministic.
        """
        def build_start_state(seq: str) -> List[Complex]:
            d_top = Domain(name="orig",
                           sequence=seq or generate_sequence(seqlen))
            s_top = Strand(name="top", domains=[d_top])
            s_bot = s_top.C
            startTop = Complex(strands=[s_top], structure=".")
            startBot = Complex(strands=[s_bot], structure=".")
            return [startTop, startBot]

        def options_factory(num_sims: int, seq: str) -> Options:
            o = self.create_config(rate_model, build_start_state(seq), seed)
            o.simulation_mode = Literals.first_step
            o.num_sims = num_sims
            o.simulation_time = 1e-6
            return o

        num_sims = 4
        seq = generate_sequence(seqlen)

        sim = MergeSim()
        sim.setNumOfThreads(num_sims)
        with raising_random:
            sim.setOptionsFactory2(options_factory, num_sims,
                                   "" if random_seq else seq)
        sim.setOptionsFactory2(options_factory, num_sims, seq)
        if wrong_mode:
            sim.setPassageMode()
        else:
            sim.setFirstStepMode()
        with raising_mode:
            sim.run()


# ==============================================================================


if __name__ == "__main__":
    toehold_seq = "GTGGGT"
    bm_design_B = "ACCGCACGTCACTCACCTCG"
    toehold_extra = "TTT"
    structure = "..(.((.....)).).....((((((+))))))((((((((((((((((((((+))))))))))))))))))))"
    seed = np.random.randint(int(1e10))
    Test_Determinism().test_trajectories(
        "DNA23Arrhenius", toehold_seq, bm_design_B, toehold_extra, structure,
        seed)
