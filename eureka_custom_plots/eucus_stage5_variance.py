"""
OVERVIEW:

WILL BE FILLED IN!

"""

"""
PSEUDOCODE TODO:
- list every filename connected to data type
- make computations for every file
- White noise only necessary for one file!
"""
import h5py as h5
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from util import rc_setup

# GLOBALS
DATA_DIR = "data"
FILE_NAME = "S5_wasp39b_ap4_bg8_Table_Save_ch00.txt"
PLOT_DIR = "plots/stage5_variance"
PLOT_TYPE = "png"


def main():
    resid = lc_data(f"{DATA_DIR}/{FILE_NAME}")
    rms, bins, stderr = allan_values(resid, binstep=1)
    allan_plot(bins, rms, stderr)


def lc_data(filename):
    """DOC!"""
    data = pd.read_csv(filename, sep=" ", comment="#")

    return data["residuals"]


def allan_values(residuals, binstep=1, maxbins=None):
    """DOC!"""
    if maxbins is None:
        maxbins = residuals.shape[0] / 10

    binsz = np.arange(1, maxbins + binstep, step=binstep, dtype=int)
    nbins = np.zeros(binsz.size, dtype=int)
    rms = np.zeros(binsz.size)

    for i in range(binsz.size):
        nbins[i] = int(np.floor(residuals.size / binsz[i]))
        bindata = np.ma.zeros(nbins[i], dtype=float)
        # bin data
        # ADDED INTEGER CONVERSION, mh 01/21/12
        for j in range(nbins[i]):
            bindata[j] = np.ma.mean(residuals[j * binsz[i]:(j + 1) * binsz[i]])
        # get rms
        rms[i] = np.sqrt(np.ma.mean(bindata ** 2))

    temp = np.sqrt(nbins / (nbins - 1.))
    stderr = np.std(residuals) / np.sqrt(binsz) * temp

    return rms, binsz, stderr


def allan_plot(bins, rms, stderr):
    """DOC"""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot wavelength-channel RMS
    ax.plot(bins, rms, ls="-", lw=2, c="black", alpha=0.5)

    # Plot white noise expected RMS
    ax.plot(bins, stderr, ls="--", lw=2.5, c="tab:red")

    ax.set(
        xscale="log", xlabel="Bin size (# of integrations)",
        yscale="log", ylabel="RMS"
    )

    plt.savefig(f"{PLOT_DIR}/comb_allan_variance.{PLOT_TYPE}", dpi=600)
    plt.show()


if __name__ == "__main__":
    rc_setup()
    main()
