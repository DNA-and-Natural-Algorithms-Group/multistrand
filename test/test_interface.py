# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from random import randrange

import pytest

from multistrand.objects import Domain, Strand, Complex, StopCondition
from multistrand.options import Options
from multistrand._options.interface import FirstStepResult
from multistrand._options.constants import OptionsConstants
from multistrand.system  import SimSystem


@pytest.mark.parametrize("kinetics",
                         ["JSMetropolis25", "DNA23Metropolis", "DNA23Arrhenius"])
@pytest.mark.parametrize("num_simulations", [100])
class Test_Interface:
    """
    This test exercises the Python interface.
    """

    def test_first_step_mode(self, kinetics: str, num_simulations: int,
                             capfd: pytest.CaptureFixture):
        print("Creating options...")
        d1 = Domain(name="d1", sequence="GTTGGTTTGTGTTTGGTGGG")
        s1 = Strand(name="s1", domains=[d1])
        c1 = Complex(name="c1", strands=[s1], structure=".")
        c2 = Complex(name="c2", strands=[s1.C], structure=".")
        c3 = Complex(name="c3", strands=[s1, s1.C], structure="(+)")

        sc_rev = StopCondition("REVERSE", [(c1, 2, 0), (c2, 2, 0)])
        sc_for = StopCondition("END", [(c3, 4, 6)])
        o = Options(simulation_mode='First Step', num_simulations=num_simulations,
                    simulation_time=0.5, start_state=[c1, c2],
                    stop_conditions=[sc_rev, sc_for])
        getattr(o, kinetics)()
        o.initial_seed = randrange(-2147483648, 2147483647)

        print("\nCreating SimSystem...")
        s = SimSystem(o)

        print("\nStarting simulation...")
        s.start()

        print("\nInterrogating interface...")
        print(o.interface)
        results = o.interface.results
        print(results)

        assert len(results) == len(o.interface.end_states) == num_simulations
        for res in results:
            assert isinstance(res, FirstStepResult)
            assert isinstance(res.seed, int)
            assert res.com_type == OptionsConstants().STOPRESULT["Normal"]
            assert res.tag in [sc.tag for sc in [sc_rev, sc_for]]
            assert res.time > 0
            assert res.collision_rate > 0

        capture = capfd.readouterr()
        assert all(capture.out.find(s) > 0 for s in
                   ["trajectories completed", "Seed", "Result",
                    "Completion Time", "Collision Rate", "Completion Tag",
                    "Normal"])
        print(capture.out)
