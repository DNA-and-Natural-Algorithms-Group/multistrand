# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

import numpy as np
from types import ModuleType
from typing import IO, Optional

import pytest

from multistrand.objects import Domain
import multistrand.utils.thermo as thermo


class Test_Dissoc_Rate:
    """
    This test compares Multistrand methods to compute dissociation rates for
    duplexes.

    We can either compute dissociation rate k- directly or compute k+ and use
    the Boltzmann equation k+ / k- = exp( - dG / RT ) to infer k-. dG is
    obtained from NUPACK. As it turns out, using either method is equivalent,
    both for Arrhenius and Metropolis models.
    """

    @staticmethod
    def comparison(seq: str, tutorials: ModuleType, f: Optional[IO] = None):
        seqC = Domain(name="top", sequence=seq).C.sequence
        temp = 25.0 + 40 + thermo.C2K
        sodium = 1.0

        hline = 60 * "="
        if f is not None:
            f.write(str(seq) + "   ")
            f.write(str("%0.3g" % (temp)) + "   ")

        from tutorials.compute.anneal import compute as computeA
        simA = computeA(seq, temp, sodium)
        predictedA = simA.results
        print(predictedA)
        print(hline)

        from tutorials.compute.dissociation import compute as computeD
        simD = computeD(seq, temp, sodium)
        predictedD = simD.results
        print(predictedD)
        print(hline)

        # dotparen = "("*len(seq) + "+" + ")"*len(seq)
        model = thermo.Model(material="DNA", kelvin=temp, ensemble="some", sodium=1.0, magnesium=0.0)
        dG = thermo.pfunc(strands=[seq, seqC], model=model)
        print(f"temp {temp}")
        print(f"dG = {dG}")

        k1_A = predictedA.k1() * np.exp(dG / (thermo.GAS_CONSTANT * temp))
        k1_D = predictedD.k1()
        log_diff = np.abs(np.log10(k1_A) - np.log10(k1_D))
        print(f"{k1_A = :.3g}, {k1_D = :.3g}, {log_diff = :.3g}")

        if f is not None:
            f.write(str("%0.3g" % np.log10(k1_A)) + "    ")
            f.write(str("%0.3g" % np.log10(predictedD.k1())) + "    ")
            f.write(str("%0.3g" % log_diff) + "    ")
            f.write("\n")
            f.flush()
        return (simA, simD), (k1_A, k1_D, log_diff)

    @pytest.mark.parametrize("seq", ["AGCTGATTC"])
    def test_comparison(self, seq: str, tutorials: ModuleType,
                        capfd: pytest.CaptureFixture):
        # run simulation scripts and extract rate estimates
        (simA, simD), log_ks = self.comparison(seq, tutorials)

        # check that simulation environments were configured identically, except
        # for options concerning the reaction direction
        optA, optD = simA.factory.new(0), simD.factory.new(0)
        optD.simulation_mode = optA.simulation_mode
        optD.simulation_time = optA.simulation_time
        optD._start_state = optA._start_state.copy()
        optD._stop_conditions = optA._stop_conditions.copy()
        assert optA == optD

        # check consistency of rate estimates
        assert all(isinstance(k, float) and k > 0 for k in log_ks)
        log_diff = log_ks[2]
        assert log_diff < .2

        # check expected script output
        capture = capfd.readouterr()
        print(capture.out)
        assert capture.out != ""
        assert capture.err == ""
        assert all(capture.out.find(s) > 0 for s in
                   ["FirstStepRate", "FirstPassageRate", "Found",
                    "successful trials", "nForward", "nReverse"])


def main():
    # For each column, compute the rate and print it to a file.
    resultFileName = "dissociation_comparison.txt"
    with open(resultFileName, 'w+') as f:
        f.write("Seq    temp     predicted-D    predicted-A   difference \n\n")

        seq = "ATCCCAATCAACACCTTTCCTA"
        Test_Dissoc_Rate.comparison(seq, f)


if __name__ == '__main__':
    main()
