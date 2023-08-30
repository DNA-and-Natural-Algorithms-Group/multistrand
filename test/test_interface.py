# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from random import randrange
from typing import List

import pytest

from multistrand.options import Options, Literals, EnergyType
from multistrand._options.interface import FirstStepResult, OptionsConstants
from multistrand.objects import (
    Domain, ComplementaryDomain, Strand, Complex, StopCondition)
from multistrand.system  import SimSystem, calculate_energy


class Test_Interface:
    """
    This test exercises the Python interface.
    """
    @staticmethod
    def small_system():
        domains = [Domain(name="d1", sequence="ACTTG")]
        domains.append(ComplementaryDomain(domains[0]))
        strands = [Strand("s1", domains=[domains[0]]),
                   Strand("s2", domains=[domains[1]])]
        complexes = [Complex(name="c1", strands=[strands[0]], structure="....."),
                     Complex(name="c2", strands=[strands[1]], structure="....."),
                     Complex(name="c3", strands=strands, structure="(((((+)))))")]
        conditions = [StopCondition("REVERSE", [(complexes[0], 2, 0),
                                                (complexes[1], 2, 0)]),
                      StopCondition("END", [(complexes[2], 4, 1)])]
        return domains, strands, complexes, conditions

    @staticmethod
    def big_system():
        d1 = Domain(name="d1", sequence="GTTGGTTTGTGTTTGGTGGG")
        s1 = Strand(name="s1", domains=[d1])
        c1 = Complex(name="c1", strands=[s1], structure=".")
        c2 = Complex(name="c2", strands=[s1.C], structure=".")
        c3 = Complex(name="c3", strands=[s1, s1.C], structure="(+)")
        sc_rev = StopCondition("REVERSE", [(c1, 2, 0), (c2, 2, 0)])
        sc_for = StopCondition("END", [(c3, 4, 6)])
        return d1, s1, c1, c2, c3, sc_rev, sc_for

    @staticmethod
    def check_duplicates(seq: List[Domain | Strand | Complex | StopCondition]):
        """
        Helper function to see if we duplicated the 'id' attribute or 'tag'
        attribute in an interable.
        """
        ob_type = seq[0].__class__.__name__
        assert len(seq) == len(set([
            f"{i.__class__.__name__}[{(i.tag if isinstance(i, StopCondition) else i.id)}]"
            for i in seq])), f"Duplicate {ob_type} instances were added"

    def test_base_objects(self):
        """
        This test creates all the basic objects.
        """
        domains, strands, complexes, conditions = self.small_system()
        self.check_duplicates(domains)
        self.check_duplicates(strands)
        self.check_duplicates(complexes)
        self.check_duplicates(conditions)

    @pytest.mark.parametrize(
        "mode", [Literals.first_passage_time, Literals.first_step])
    @pytest.mark.parametrize(
        "kinetics", ["JSMetropolis25", "DNA23Metropolis", "DNA23Arrhenius"])
    @pytest.mark.parametrize("num_init", [5])
    @pytest.mark.parametrize("num_run", [5])
    def test_options_simsystem(self, kinetics: str, mode: Literals,
                               num_init: int, num_run: int):
        """
        Create an `Options` object using some common values and simple
        start/stop states, and attempt its reuse across multiple `SimSystem`
        objects.
        """
        _, _, complexes, conditions = self.small_system()

        o = Options()
        o.simulation_mode = mode
        o.substrate_type  = Literals.substrateDNA
        getattr(o, kinetics)()
        o.num_simulations = 3
        o.output_interval = 1
        o.start_state = complexes
        if mode == Literals.first_passage_time:
            o.stop_conditions = conditions
            o.simulation_time = 0.5
        elif mode == Literals.first_step:
            o.simulation_time = 1e-4
        else:
            raise ValueError

        for _ in range(num_init):
            s = SimSystem(o)
        for _ in range(num_run):
            assert isinstance(o, Options)
            assert isinstance(s, SimSystem)
            print("Start simulation...")
            s.start()
            print("Interface:")
            print(o.interface)
            print("Energy calculation:",
                  calculate_energy([complexes[1]], o, EnergyType.complex))
            print()

    @pytest.mark.parametrize(
        "kinetics", ["JSMetropolis25", "DNA23Metropolis", "DNA23Arrhenius"])
    @pytest.mark.parametrize("num_simulations", [100])
    def test_first_step_mode(self, kinetics: str, num_simulations: int,
                             capfd: pytest.CaptureFixture):
        print("Creating options...")
        _, _, c1, c2, _, sc_rev, sc_for = self.big_system()
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
            assert res.com_type == OptionsConstants.STOPRESULT["Normal"]
            assert res.tag in [sc.tag for sc in [sc_rev, sc_for]]
            assert res.time > 0
            assert res.collision_rate > 0

        capture = capfd.readouterr()
        assert all(capture.out.find(s) > 0 for s in
                   ["trajectories completed", "Seed", "Result",
                    "Completion Time", "Collision Rate", "Completion Tag",
                    "Normal"])
        print(capture.out)
