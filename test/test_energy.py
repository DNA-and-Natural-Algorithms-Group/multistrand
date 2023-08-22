# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from collections import defaultdict
from functools import partial
from multiprocess import Pool, cpu_count
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Optional

import numpy as np
import pytest

from multistrand.options import Options, Energy_Type
from multistrand.system import energy
from multistrand.objects import Complex, Strand
import multistrand.utils.thermo as thermo


def split_tasks(n_tasks: int, n_workers: int) -> List[np.ndarray]:
    """
    Partition the index set of `n_tasks` into at most `n_workers` blocks.
    """
    n_rem = n_tasks % (n_workers - 1)
    tasks_div = np.arange(n_tasks - n_rem)
    tasks_rem = np.arange(n_tasks - n_rem, n_tasks)
    blocks = np.split(tasks_div, n_workers - 1) + [tasks_rem]
    assert len(blocks) <= n_workers
    return blocks


class Test_SingleStrandEnergy:
    """
    Compare the thermodynamic scores between Multistrand and Nupack, for a set
    of single-stranded complexes.
    """
    # subsampling of input file, defined per category
    examples_fraction: float = 1.0  # all examples when 1.0
    examples_min: int = 1000

    # runtime config
    num_workers = max(2, cpu_count() - 2)

    @pytest.mark.parametrize("rel_tol", [1e-6])
    @pytest.mark.parametrize(
        "examples_file", [Path(__file__).parent / 'testSetSS.txt'])
    def test_energy(self, examples_file: Path, rel_tol: float):
        complexes = self.load_complexes(examples_file)
        with Pool() as pool:
            for category, (seqs, structs) in complexes.items():
                print(f"{category}: {len(seqs)}")
                blocks = [(seqs[ix], structs[ix])
                          for ix in split_tasks(len(seqs), self.num_workers)]
                pool.map(partial(self.compare_energies, rel_tol, category),
                         blocks)

    def load_complexes(self, path: Path) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        dataset = defaultdict(list)
        category: Optional[str] = None
        # read complexes from file
        with open(path) as f:
            for l in f:
                if l.startswith('>'):
                    category = l[1:].strip()
                else:
                    assert category is not None
                    dataset[category].append(l.strip())
        # parse and subsample
        complexes = {}
        for category, samples in dataset.items():
            seqs, structs = np.array(samples[0::2]), np.array(samples[1::2])
            N = len(seqs)
            n = N if N <= self.examples_min else int(self.examples_fraction * N)
            idx = np.arange(N)
            np.random.shuffle(idx)
            idx = idx[:n]
            complexes[category] = (seqs[idx], structs[idx])
        return complexes

    @staticmethod
    def create_config() -> Options:
        opt = Options(verbosity=0)
        opt.DNA23Metropolis()
        return opt

    @classmethod
    def compare_energies(cls, rel_tol: float, category: str,
                         complexes: Tuple[Iterable[str], Iterable[str]]) -> None:
        opt = cls.create_config()
        model = thermo.Model(opt)
        i = 0
        for seq, struct in zip(*complexes):
            assert len(seq) == len(struct)
            e_nupack = thermo.energy([seq], struct, model=model)
            c_multistrand = Complex(
                strands=[Strand(name="hairpin", sequence=seq)], structure=struct)
            e_multistrand = energy(
                [c_multistrand], opt, Energy_Type.Complex_energy)
            assert np.allclose(e_nupack, e_multistrand, rtol=rel_tol), \
                f"category = {category}, seq = {seq}, struct = {struct}"


if __name__ == "__main__":
    test = Test_SingleStrandEnergy()
    test.examples_fraction = 0.01
    test.test_energy(Path(__file__).parent / 'testSetSS.txt', rel_tol=1e-6)
