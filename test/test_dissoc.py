# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

import math
import numpy as np
from types import ModuleType
from typing import IO, Optional

import pytest

import multistrand
from multistrand.options import Options
from multistrand.objects import Domain
from multistrand.utils import GAS_CONSTANT, C2K
from nupack import pfunc


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
        from tutorials.compute.anneal import compute as computeAnneal
        from tutorials.compute.dissociation import compute as computeDissociation

        seqC = Domain(name="top", sequence=seq).C.sequence
        temp = 25.0 + 40 + C2K
        sodium = 1.0

        hline = 60 * "="
        if f is not None:
            f.write(str(seq) + "   ")
            f.write(str("%0.3g" % (temp)) + "   ")

        simA = computeAnneal(seq, temp, sodium)
        predictedA = simA.results
        print(predictedA)
        print(hline)

        simD = computeDissociation(seq, temp, sodium)
        predictedD = simD.results
        print(predictedD)
        print(hline)

        # dotparen = "("*len(seq) + "+" + ")"*len(seq)
        dG = pfunc([seq, seqC], [1, 2], T=temp - C2K, material="dna")
        print()
        print(f"{temp = }")
        print(f"{dG = }")

        k1_A = predictedA.k1() * math.exp(dG / (GAS_CONSTANT * temp))
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
        optA, optD = simA.factory.new(0), simD.factory.new(0)
        log_diff = log_ks[2]

        # ensure that simulations were configured identically
        assert optA.rate_method == optD.rate_method
        assert optA.unimolecular_scaling == optD.unimolecular_scaling
        assert optA.bimolecular_scaling == optD.bimolecular_scaling
        assert (optA.lnAStack, optA.EStackEnd) == (optD.lnAStack, optD.EStackEnd)

        # check consistency of rate estimates
        assert all(isinstance(k, float) and k > 0 for k in log_ks)
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
