"""
Read T-Walk HDF5 results
========================
"""

import sys

import arviz as az
import h5py
import numpy as np
import matplotlib.pyplot as plt
import corner


class TWalk(h5py.File):
    """
    File object for T-Walk results.
    """
    def __init__(self, *args, **kwargs):
        """
        Initialize T-Walk specific attributes
        """
        super().__init__(*args, **kwargs)
        self.data = self[list(self.keys())[-1]]
        self.valid = np.array(self.data[f"chain_isvalid"][:])
        self.chain = self.get_array("chain")
        self.n_chains = max(self.chain) + 1
        self.mults = self.get_chains("mult")

    def get_array(self, key):
        """
        :returns: Array of valid data for a key
        """
        return np.array(self.data[key][:])[self.valid == 1]

    def get_chains(self, key):
        """
        :returns: List of arrays of valid data for a key, each one from a single walker
        """
        arr = self.get_array(key)
        return [arr[self.chain == i] for i in range(self.n_chains)]

    def get_equally_weighted_chains(self, key):
        """
        :returns: List of arrays of valid data for a key, each one from a single walker
        """
        chains = self.get_chains(key)
        return [equally_weight(c, m) for c, m in zip(chains, self.mults)]

    def get_param_names(self):
        """
        :returns: Guess names of paramaters of interest in chains
        """
        return [k for k in self.data.keys()
                if "::primary_parameters::" in k and not k.endswith("_isvalid")]

    def to_arviz(self):
        """
        :returns: Arviz data object for chain
        """
        dict_ = {k.split("::")[-1]: equal_length_vstack(
            self.get_equally_weighted_chains(k)) for k in self.get_param_names()}
        return az.dict_to_dataset(dict_)


def equally_weight(chain, mult):
    """
    :returns: Repeat states according to their multiplicity
    """
    e = []
    for c, m in zip(chain, mult):
        e += [c] * m
    return e


def equal_length_vstack(chains):
    """
    :returns: Stack unequal length chains by making them all no longer than the shortest
    """
    n = min([len(c) for c in chains])
    return np.vstack([c[:n] for c in chains])


if __name__ == "__main__":

    h5_name = sys.argv[1]

    with TWalk(h5_name, "r") as f:
        data = f.to_arviz()

    az.plot_trace(data, compact=True)
    plt.savefig("trace.pdf")

    az.plot_pair(data, kind='kde')
    plt.savefig("pair.pdf")

    az.plot_posterior(data)
    plt.savefig("posterior.pdf")

    az.plot_autocorr(data)
    plt.savefig("autocorr.pdf")

    corner.corner(data)
    plt.savefig("corner.pdf")

    print(az.summary(data))
    print(az.rhat(data))
