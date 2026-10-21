import warnings
import numpy as np

from ase import Atoms
from dscribe.descriptors import MBTR


def fp_mbtr(atoms: Atoms, periodic: bool, normalization: str) -> list:
    mbtr = MBTR(
        species=list(set(atoms.get_atomic_numbers())),
        k2={
            "geometry": {"function": "distance"},
            "grid": {"min": 0.0, "max": 5.0, "sigma": 0.1, "n": 100},
            "weighting": {"function": "inverse_square", "r_cut": 4.0, "scale": 0.5, "threshold": 1e-3},
        },
        k3={
            "geometry": {"function": "cosine"},
            "grid": {"min": -1.0, "max": 1.0, "sigma": 0.1, "n": 50},
            "weighting": {"function": "exp", "r_cut": 3.0, "scale": 0.3, "threshold": 1e-3},
        },
        periodic=periodic,
        sparse=False,
        flatten=False,
        normalization=normalization
    )
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", message=".*invalid value encountered in true_divide.*",
                                category=RuntimeWarning)
        mbtr_fp = mbtr.create(system=atoms, n_jobs=1, only_physical_cores=False)
    k2 = np.sum(mbtr_fp["k2"], axis=(0, 1))
    k3 = np.sum(mbtr_fp["k3"], axis=(0, 1, 2))
    # 把k1, k2, k3中所有的Nan值赋为0
    k2[np.isnan(k2)] = 0
    k3[np.isnan(k3)] = 0
    return np.concatenate((k2, k3)).tolist()