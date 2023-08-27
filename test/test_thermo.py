# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from itertools import repeat

import numpy as np
import nupack

from multistrand.utils.thermo import (
    energy, defect, prob, sample, mfe, pairs,
    complex_free_energy, meltingTemperature, C2K)


def test_thermo_wrapper():
    """
    Usage examples for the `utils.thermo` wrapper around NUPACK utilities [1].

    NOTE: "In NUPACK 4, pseudoknots are excluded from the structural ensemble."

    [1] https://docs.nupack.org/
    """
    rna_seq = ['GGGCUGUUUUUCUCGCUGACUUUCAGCCCCAAACAAAAAAUGUCAGCA']
    dna_pair = ['ACGT','GCTT']
    dna_seqs = ['AGTCTAGGATTCGGCGTGGGTTAA',
                'TTAACCCACGCCGAATCCTAGACTCAAAGTAGTCTAGGATTCGGCGTG',
                'AGTCTAGGATTCGGCGTGGGTTAACACGCCGAATCCTAGACTACTTTG']
    dna_struct = '+'.join(['((((((((((((((((((((((((',
                           '))))))))))))))))))))))))((((((((((((.(((((((((((',
                           '........................))))))))))).))))))))))))'])

    dG_ensemble = complex_free_energy(dna_seqs, material='dna')
    assert np.allclose(dG_ensemble, -62.6645)
    print(f"Complex free energy of 3 DNA strands: {dG_ensemble:.3f} kcal/mol")

    dG_ensemble = complex_free_energy(rna_seq, material='rna')
    assert np.allclose(dG_ensemble, -11.8452)
    print(f"Complex free energy of 1 RNA strand: {dG_ensemble:.3f} kcal/mol")

    dG = energy(dna_seqs, dna_struct, material='dna')
    assert np.allclose(dG, -61.7927)
    print(f"Free energy of a DNA structure: {dG:.3f} kcal/mol")

    n = defect(dna_struct, dna_seqs, material='dna')
    assert np.allclose(n, 0.06914)
    print(f"Normalised complex ensemble defect for a DNA structure: {n:.3%}")

    p = prob(dna_seqs, dna_struct, material='dna')
    assert np.allclose(p, 7.992e-05)
    print(f"Probability of a DNA structure: {p:%}")

    structs = sample(dna_seqs, 3, material='dna')
    print("Boltzmann sample of DNA structures:")
    for s in structs:
        assert isinstance(s, nupack.Structure)
        assert ((np.array(s.dp()) == '+') == (np.array(dna_struct) == '+')).all()
        print(f"  {s}")

    s0 = mfe(dna_seqs, material='dna')[0]
    assert isinstance(s0, nupack.thermo.StructureEnergy)
    print("MFE structure:")
    print(f"  energy = {s0.energy:.3f} kcal/mol")
    print(f"  {s0.structure}")

    T = meltingTemperature(dna_seqs[0])
    print(f"Melting temperature of a DNA duplex: {T - C2K:.3f} C")

    pm = pairs(dna_pair, material='dna')
    assert isinstance(pm, nupack.PairMatrix)
    pmA = pm.to_array()
    assert pmA.shape == tuple(repeat(sum(map(len, dna_pair)), 2))
    assert (0 <= pmA).all() and (pmA <= 1).all()
    print("Pair probabilities of DNA strands:")
    print("  " + str(pm).replace("\n", "\n  "))
