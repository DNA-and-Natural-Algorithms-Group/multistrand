"""
Testing loop internals.
"""

from typing import List, Optional
from functools import reduce

import pytest

from multistrand.objects import Strand, Complex, Domain, StopCondition
from multistrand.options import Options, Literals
from multistrand.system import SimSystem
from multistrand.experiment import makeComplex


# ==============================================================================


@pytest.mark.parametrize("kinetics",
                         ["JSMetropolis25", "DNA23Metropolis", "DNA23Arrhenius"])
class Test_InitialInfo:
    """
    Check the successful and deterministic construction of the initial SimSystem
    state for a variety of system specifications.
    """
    @staticmethod
    def create_options(start_state: List[Complex], stop_complex: Complex,
                       simMode=Literals.first_passage_time):
        return Options(
            simulation_mode=simMode,
            parameter_type="Nupack", substrate_type="DNA",
            temperature=273.15 + 25.0,
            num_simulations=10, simulation_time=0.00001,
            # rate_scaling='Calibrated',
            join_concentration=1.0, rate_method="Metropolis",
            start_state=start_state,
            stop_conditions=[StopCondition(
                "CLOSED", [(stop_complex, Literals.dissoc_macrostate, 2)])])

    @staticmethod
    def print_initialInfo(opt: Options):
        opt.verbosity = 1
        sim = SimSystem(opt)
        sim.initialInfo()

    @classmethod
    def compare_initialInfo(cls, opt: Options, kinetics: str,
                            capfd: Optional[pytest.CaptureFixture]):
        """
        Test helper. Call with `capfd=None` to debug without Pytest.
        """
        # set kinetic parameters
        getattr(opt, kinetics)()

        # initialise full simulator state twice
        _ = None if capfd is None else capfd.readouterr()
        cls.print_initialInfo(opt)
        capture1 = None if capfd is None else capfd.readouterr()
        if capfd is None:
            print(50 * '-')
        cls.print_initialInfo(opt)
        capture2 = None if capfd is None else capfd.readouterr()

        # compare & check outputs
        if capfd is not None:
            print(capture1.out)
            print(50 * '-')
            assert capture1.out == capture2.out != ""
            assert capture1.err == capture2.err == ""
            assert capture1.out.startswith("Complex ")
            assert all(capture1.out.find(s) > 0 for s in
                       ["seq", "struc", "energy-ms", "energy-nu", "***",
                        "JoinFlux", "joinrate"])

    @classmethod
    def debug_single_test(cls, test_id: str, *args):
        kinetics = "JSMetropolis25"
        getattr(cls(), f"test_{test_id}")(kinetics, None, *args)

    # --------------------------------------------------------------

    @pytest.mark.parametrize("toehold_seq, domain_seq",
                             [("CCCC", "CATTAAC"), ("CC", "CAAC")])
    @pytest.mark.parametrize(
        "start_struct, stop_struct",
        [("((+))", "..+.."),
         pytest.param("..+..", "((+))",
                      marks=pytest.mark.skip(
                          reason="Unconnected strand in initialized complex. "
                          "| Segmentation fault in OpenLoop::getOpenInfo()"))])
    def test_0(self, kinetics: str, capfd: pytest.CaptureFixture,
               toehold_seq: str, domain_seq: str,
               start_struct: str, stop_struct: str):
        """
        A fully hybridized strand.
        """
        toehold = Domain(name="toehold", sequence=toehold_seq,
                            length=len(toehold_seq))
        branch_migration = Domain(name="branch_migration", sequence=domain_seq,
                                    seq_length=len(domain_seq))
        incoming = branch_migration.C + toehold.C
        substrate = toehold + branch_migration
        start_complex = Complex(strands=[incoming, substrate], structure=start_struct)
        stop_complex = Complex(strands=[incoming, substrate], structure=stop_struct)
        opt = self.create_options([start_complex], stop_complex)
        self.compare_initialInfo(opt, kinetics, capfd)

    def test_0C(self, kinetics: str, capfd: pytest.CaptureFixture):
        domain_seq, domain_seq2 = "AGT", "GTA"
        left = Domain(name="branch_migration", sequence=domain_seq, seq_length=3)
        right = Domain(name="branch_migration2", sequence=domain_seq2, seq_length=3)
        incoming = left + right
        substrate = incoming.C
        start_complex1 = Complex(strands=[incoming], structure="..")
        start_complex2 = Complex(strands=[substrate], structure="..")
        stop_complex = Complex(strands=[incoming, substrate], structure="((+))")
        opt = self.create_options([start_complex1, start_complex2], stop_complex)
        # FIXME: Obsolete options?
        # opt.rate_scaling = 'Calibrated'
        opt.join_concentration = 1e-9
        self.compare_initialInfo(opt, kinetics, capfd)

    @pytest.mark.parametrize("top0, top1, toehold, bottom",
                             [("ACT", "GAC", "TG", ""),
                              ("ACT", "GAC", "TG", "TG")])
    def test_1(self, kinetics: str, capfd: pytest.CaptureFixture,
               top0: str, top1: str, toehold: str, bottom: str):
        """
        Open loop, w/ and w/o initialization penalty.
        """
        right_d = Domain(name="toehold0", sequence=top0, length=len(top0))
        left_d = Domain(name="toehold1", sequence=top1, length=len(top1))
        branch0 = Domain(name="branch_migration", sequence=bottom, seq_length=len(bottom))
        toehold = Domain(name="oToehold", sequence=toehold, length=len(toehold))
        substrate = toehold.C + branch0.C + toehold.C
        left = toehold + left_d
        right = right_d + toehold
        # Note that "+" is used to indicate strand breaks.
        # So the initial structures represent the incoming strand bound by its toehold,
        # and we'll see that either it completes strand displacement, or it dissociates.
        start_complex = Complex(strands=[left, right, substrate], structure="(.+.(+).)")
        stop_complex = Complex(strands=[left, right, substrate], structure="..+..+...")
        opt = self.create_options([start_complex], stop_complex)
        self.compare_initialInfo(opt, kinetics, capfd)

    @pytest.mark.parametrize("top, bottom",
                             [("TTT", "AAA"), ("TGG", "CCA")])
    @pytest.mark.parametrize("start_struct",
                             ["(..+..)", ".(.+.).", "..(+).."])
    def test_2(self, kinetics: str, capfd: pytest.CaptureFixture,
               top: str, bottom: str, start_struct: str):
        """
        Simple test.
        """
        seq = top + bottom
        strands = [Domain(name=f"toehold{i}", sequence=seq[i], length=1)
                   for i in range(len(top) + len(bottom))]
        substrate = reduce(Domain.__add__, strands[3:])
        invading = reduce(Domain.__add__,strands[:3])
        start = Complex(strands=[substrate, invading], structure=start_struct)
        stop = Complex(strands=[substrate, invading], structure="...+...")
        opt = self.create_options([start], stop)
        self.compare_initialInfo(opt, kinetics, capfd)

    def test_3(self, kinetics: str, capfd: pytest.CaptureFixture):
        """
        Bimolecular step.
        """
        strand_seq = "CTGA"
        num_traj = 10
        onedomain = Domain(name="itall", sequence=strand_seq)
        top = Strand(name="top", domains=[onedomain])
        bot = top.C
        # Note that the structure is specified to be single stranded, but this
        # will be over-ridden when Boltzmann sampling is turned on.
        start_complex_top = Complex(strands=[top], structure=".")
        start_complex_bot = Complex(strands=[bot], structure=".")
        start_complex_top.boltzmann_count = num_traj
        start_complex_bot.boltzmann_count = num_traj
        start_complex_top.boltzmann_sample = True
        start_complex_bot.boltzmann_sample = True
        # Turns Boltzmann sampling on for this complex and also does sampling
        # more efficiently by sampling 'num_traj' states.

        # Stop when the exact full duplex is achieved. (No breathing!)
        success_complex = Complex(strands=[top, bot], structure="(+)")
        success_stop_condition = StopCondition(
            "SUCCESS", [(success_complex, Literals.exact_macrostate, 0)])
        # Declare the simulation unproductive if the strands become
        # single-stranded again.
        failed_complex = Complex(strands=[top], structure=".")
        failed_stop_condition = StopCondition(
            "FAILURE", [(failed_complex, Literals.dissoc_macrostate, 0)])
        opt = Options(simulation_mode="First Step",
                      parameter_type="Nupack", substrate_type="DNA",
                      rate_method="Metropolis",
                      num_simulations=num_traj, simulation_time=1.0,
                      dangles="Some", temperature=273.15 + 25.0,
                      # rate_scaling="Calibrated",
                      useArrRates=True, verbosity=0)
        opt.start_state = [start_complex_top, start_complex_bot]
        opt.stop_conditions = [success_stop_condition, failed_stop_condition]
        self.compare_initialInfo(opt, kinetics, capfd)

    @pytest.mark.skip(reason="Failed assertion in StrandOrdering::addOpenLoop()")
    def test_4(self, kinetics: str, capfd: pytest.CaptureFixture):
        """
        Lightbulb.
        """
        toehold_t = "CTGC"
        toehold_dd = "CATATC"
        domain_R = "CATTAAC"
        toehold = Domain(name="toehold", sequence=toehold_t, length=6)
        toehold_2 = Domain(name="toehold", sequence=toehold_dd, length=6)
        branch_migration = Domain(name="branch_migration", sequence=domain_R, seq_length=7)
        incoming = branch_migration.C + toehold.C
        substrate = toehold + branch_migration
        # Note that "+" is used to indicate strand breaks.
        # So the initial structures represent the incoming strand bound by its toehold,
        # and we'll see that either it completes strand displacement, or it dissociates.
        start_complex = Complex(strands=[incoming, substrate], structure=".(+).")
        stop_complex = Complex(strands=[incoming, substrate], structure="..+..")
        opt = self.create_options([start_complex], stop_complex)
        self.compare_initialInfo(opt, kinetics, capfd)

    def test_5(self, kinetics: str, capfd: pytest.CaptureFixture):
        """
        Multi-loop.
        """
        top0 = "ACT"
        top1 = "GAC"
        toehold = "TG"
        bottom = "TG"
        connected = "ATA"
        right_d = Domain(name="toehold0", sequence=top0, length=3)
        left_d = Domain(name="toehold1", sequence=top1, length=3)
        branch0 = Domain(name="branch_migration", sequence=bottom, seq_length=2)
        toehold = Domain(name="oToehold", sequence=toehold, length=2)
        connect = Domain(name="cToehold", sequence=connected, length=3)
        substrate = toehold.C + branch0.C + toehold.C
        left = toehold + left_d + connect
        right = connect.C + right_d + toehold
        # Note that "+" is used to indicate strand breaks.
        # So the initial structures represent the incoming strand bound by its toehold,
        # and we'll see that either it completes strand displacement, or it dissociates.
        start_complex = Complex(strands=[left, right, substrate], structure="(.(+).(+).)")
        stop_complex = Complex(strands=[left, right, substrate], structure="...+...+...")
        opt = self.create_options([start_complex], stop_complex)
        self.compare_initialInfo(opt, kinetics, capfd)

    @pytest.mark.parametrize("toehold_seq, toehold_seq2, domain_seq",
                             [("CCC", "TTT", "AA"),
                              ("CTGC", "CAT", "CATGCTAAC")])
    def test_6(self, kinetics: str, capfd: pytest.CaptureFixture,
               toehold_seq: str, toehold_seq2: str, domain_seq: str):
        """
        Interior loop.
        """
        toehold = Domain(name="toehold", sequence=toehold_seq,
                         length=len(toehold_seq))
        toehold2 = Domain(name="toehold", sequence=toehold_seq2,
                          length=len(toehold_seq2))
        branch_migration = Domain(name="branch_migration", sequence=domain_seq,
                                  seq_length=len(domain_seq))
        incoming = toehold2.C + branch_migration.C + toehold.C
        substrate = toehold + branch_migration + toehold2
        start_complex = Complex(strands=[incoming, substrate], structure="(.(+).)")
        stop_complex = Complex(strands=[incoming, substrate], structure="...+...")
        opt = self.create_options([start_complex], stop_complex)
        self.compare_initialInfo(opt, kinetics, capfd)

    @pytest.mark.skip(reason="Mismatched ( in input. "
                             "| Segmentation fault in Loop::firstGen()")
    def test_7(self, kinetics: str, capfd: pytest.CaptureFixture):
        """
        Bulge loop.
        """
        toehold_seq = "CTGC"
        toehold_seq2 = "CAT"
        domain_seq = "CATGCTAAC"
        toehold = Domain(name="toehold", sequence=toehold_seq, length=6)
        toehold2 = Domain(name="toehold", sequence=toehold_seq2, length=6)
        branch_migration = Domain(name="branch_migration", sequence=domain_seq, seq_length=25)
        dangle = Domain(name="branch_migration", sequence="T", seq_length=1)
        incoming = toehold2.C + branch_migration.C + toehold.C
        substrate = toehold + toehold2
        start_complex = Complex(strands=[incoming, substrate], structure="(.(+))")
        stop_complex = Complex(strands=[incoming, substrate], structure="...+..")
        opt = self.create_options([start_complex], stop_complex)
        self.compare_initialInfo(opt, kinetics, capfd)

    @pytest.mark.skip(reason="Mismatched ( in input. "
                             "| Segmentation fault in Loop::firstGen()")
    def test_8(self, kinetics: str, capfd: pytest.CaptureFixture):
        """
        Hairpin loop.
        """
        toehold_seq = "CTGC"
        domain_seq = "CATGCTACAG"
        toehold = Domain(name="toehold", sequence=toehold_seq, length=6)
        branch_migration = Domain(name="branch_migration", sequence=domain_seq, seq_length=25)
        dangle = Domain(name="branch_migration", sequence="T", seq_length=1)
        incoming = toehold + branch_migration.C + toehold.C
        start_complex = Complex(strands=[incoming], structure="(.)")
        stop_complex = Complex(strands=[incoming], structure="...")
        opt = self.create_options([start_complex], stop_complex)
        self.compare_initialInfo(opt, kinetics, capfd)

    @pytest.mark.skip(reason="Unconnected strand in initialized complex."
                             "| Segmentation fault in Loop::firstGen()")
    def test_9(self, kinetics: str, capfd: pytest.CaptureFixture):
        """
        Small open loop.
        """
        seq0 = "GTGT"
        seq1 = "T"
        branch = Domain(name="toehold0", sequence=seq0, length=3)
        toehold = Domain(name="toehold1", sequence=seq1, length=1)
        ghost = Domain(name="toeholdG", sequence="T", length=1)
        substrate = toehold + branch
        left = toehold.C + ghost
        right = branch.C + ghost
        # Note that "+" is used to indicate strand breaks.
        # So the initial structures represent the incoming strand bound by its toehold,
        # and we'll see that either it completes strand displacement, or it dissociates.
        start_complex = Complex(strands=[right, left, substrate], structure="(.+(.+))")
        stop_complex = Complex(strands=[left, right, substrate], structure="..+..+..")
        opt = self.create_options([start_complex], stop_complex)
        self.compare_initialInfo(opt, kinetics, capfd)

    def test_10(self, kinetics: str, capfd: pytest.CaptureFixture):
        """
        Half open duplex.
        """
        seq0 = "GTCACTGCTTTT"
        seq1 = "GCAGTGAC"
        dotparen1 = "..((((((....+)))))).."
        seq2 = "GTCACTGC"
        dotparen2 = "........"
        complex1 = makeComplex([seq0, seq1], dotparen1)
        complex2 = makeComplex([seq2], dotparen2)
        opt = self.create_options([complex1], complex2)
        self.compare_initialInfo(opt, kinetics, capfd)

    def test_11(self, kinetics: str, capfd: pytest.CaptureFixture):
        """
        Oscillator gate.
        """
        seq0 = "GTAAAGACCAGTGGTGTGAAGATAGGAAAGGTGTTGATTGGGATTAGGAAACC"
        seq1 = "CATCACTATCAATCATACATGGTTTCCTAATCCCAATCAACACC"
        seq2 = "CATCACTATCAATCATACATGGTTTCCTATCTTCACACCACTGG"
        struc1 = ".......((((((((((((((((((((((((((((((((((((((((((((((+....................))))))))))))))))))))))))+......................))))))))))))))))))))))"
        complex1 = makeComplex([seq0, seq1, seq2], struc1)
        complex2 = makeComplex([seq0, seq1, seq2], struc1)
        opt = self.create_options([complex1], complex2)
        self.compare_initialInfo(opt, kinetics, capfd)

    def test_12(self, kinetics: str, capfd: pytest.CaptureFixture):
        """
        Displacement situation.

        This tests a certain secondary structure for wrong Arrhenius rates.

        (...........+).....((((((+))))))  GTAATGTCGGCG+CGCCGACATTAC+GTAATG      t=0.094573 ms,  dG=-4.68 kcal/mol, ArrType= 25
        (....(....).+).....((((((+))))))  GTAATGTCGGCG+CGCCGACATTAC+GTAATG      t=0.095125 ms,  dG=-3.27 kcal/mol, ArrType= 25
        Found error
        (...........+).....((((((+))))))  GTAATGTCGGCG+CGCCGACATTAC+GTAATG      t=0.095131 ms,  dG=-4.68 kcal/mol, ArrType= 0
        (.(...).....+).....((((((+))))))  GTAATGTCGGCG+CGCCGACATTAC+GTAATG      t=0.095381 ms,  dG=-1.51 kcal/mol, ArrType= 25
        (...........+).....((((((+))))))  GTAATGTCGGCG+CGCCGACATTAC+GTAATG      t=0.095417 ms,  dG=-4.68 kcal/mol, ArrType= 25
        """
        seq0 = "GTAATGTCGGCG"
        seq1 = "CGCCGACATTAC"
        seq2 = "GTAATG"
        struc1 = "(...........+).....((((((+))))))"
        complex1 = makeComplex([seq0, seq1, seq2], struc1)
        complex2 = makeComplex([seq0, seq1, seq2], struc1)
        opt = self.create_options([complex1], complex2)
        self.compare_initialInfo(opt, kinetics, capfd)

    def test_13(self, kinetics: str, capfd: pytest.CaptureFixture):
        """
        Dissociation.
        """
        seq0 = "ACTGACTGACTG"
        seq1 = "ACTG"
        seq2 = "CATTCAGTACAGT"
        struc1 = "((((((.(....+...(+)...).)).))))"
        complex1 = makeComplex([seq0, seq1, seq2], struc1)
        opt = self.create_options([complex1], complex1)
        opt.join_concentration = 1.0e-9
        self.compare_initialInfo(opt, kinetics, capfd)


# ==============================================================================


if __name__ == "__main__":
    test_id = "0"
    test_args = ["CC", "CAAC", "..+..", "((+))"]
    Test_InitialInfo.debug_single_test(test_id, *test_args)
