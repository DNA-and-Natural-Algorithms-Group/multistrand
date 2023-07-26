
import math
import numpy as np
from types import ModuleType
from typing import IO, Optional

import pytest

import multistrand
from multistrand.objects import Domain
from multistrand.utils import GAS_CONSTANT
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
        temp = 25.0 + 273.15 + 40
        sodium = 1.0

        if f is not None:
            f.write(str(seq) + "   ")
            f.write(str("%0.3g" % (temp)) + "   ")

        predictedA = computeAnneal(seq, temp, sodium)
        print(predictedA)
        print(60 * "=")
        predictedD = computeDissociation(seq, temp, sodium)
        print(predictedD)

        # dotparen = "("*len(seq) + "+" + ")"*len(seq)
        dG = pfunc([seq, seqC], [1, 2], T=(temp - 273.15), material="dna")
        print(f"dG = {dG}")

        k1_A = predictedA.k1() * math.exp(dG / (GAS_CONSTANT * temp))
        k1_D = predictedD.k1()
        log_diff = np.abs(np.log10(k1_A) - np.log10(k1_D))
        print(f"k1_A = {k1_A:.3g}, k1_D = {k1_D:.3g}, logdiff = {log_diff:.3g}")

        if f is not None:
            f.write(str("%0.3g" % np.log10(k1_A)) + "    ")
            f.write(str("%0.3g" % np.log10(predictedD.k1())) + "    ")
            f.write(str("%0.3g" % log_diff) + "    ")
            f.write("\n")
            f.flush()
        return k1_A, k1_D, log_diff

    @pytest.mark.parametrize("seq", ["AGCTGATTC"])
    def test_comparison(self, seq: str, tutorials: ModuleType,
                        capfd: pytest.CaptureFixture):
        _, _, log_diff = log_ks = self.comparison(seq, tutorials)
        assert all(isinstance(k, float) for k in log_ks)
        assert all(k > 0 for k in log_ks)
        assert log_diff <= .9

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
