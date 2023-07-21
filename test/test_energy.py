
from collections import defaultdict
from functools import partial
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Optional

import numpy as np
import pytest

from multistrand.options import Options, Energy_Type
from multistrand.system import initialize_energy_model, energy
from multistrand.objects import Complex, Strand
import nupack


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
    @pytest.mark.parametrize("examples_file", [Path(__file__).parent / 'testSetSS.txt'])
    def test_energy(cls, examples_file: Path, rel_tol: float):
        complexes = cls.load_complexes(examples_file)
        opt = cls.create_config()
        pool = Pool()
        for category, (seqs, structs) in complexes.items():
            print(f"{category}: {len(seqs)}")
            blocks = [(seqs[ix], structs[ix])
                      for ix in split_tasks(len(seqs), cls.num_workers)]
            pool.map(partial(cls.compare_energies, opt, rel_tol, category),
                     blocks)

    @classmethod
    def load_complexes(cls, path: Path) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
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
            n = N if N <= cls.examples_min else int(cls.examples_fraction * N)
            idx = np.arange(N)
            np.random.shuffle(idx)
            idx = idx[:n]
            complexes[category] = (seqs[idx], structs[idx])
        return complexes

    @staticmethod
    def create_config() -> Options:
        opt = Options()
        opt.verbosity = 0
        opt.DNA23Metropolis()
        initialize_energy_model(opt)
        return opt

    @staticmethod
    def compare_energies(opt: Options, rel_tol: float, category: str,
                         complexes: Tuple[Iterable[str], Iterable[str]]) -> None:
        for seq, struct in zip(*complexes):
            assert len(seq) == len(struct)
            e_nupack = nupack.energy([seq], struct, material='dna')
            c_multistrand = Complex(
                strands=[Strand(name="hairpin", sequence=seq)], structure=struct)
            e_multistrand = energy(
                [c_multistrand], opt, Energy_Type.Complex_energy)
            assert np.allclose(e_nupack, e_multistrand, rtol=rel_tol), \
                f"category = {category}, seq = {seq}, struct = {struct}"
