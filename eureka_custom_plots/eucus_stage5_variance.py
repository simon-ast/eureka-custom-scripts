"""
OVERVIEW:

Customised plotting routine for Stage 5 variance of lightcurve fit. This
was also presented in Lustig-Yaeger et al. (2022), Fig. 5. Shows the
possibility of correlated noise in the form of Allan-Variance-Plots.

"""

"""
PSEUDOCODE TODO:
- Insert a logging solution that can be turned on and off
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from util import rc_setup

# GLOBALS
DATA_DIR = "data/stage5_lcfit_results"
PLOT_DIR = "plots/stage5_variance"
PLOT_TYPE = "png"


def main():
    # Specify a list of fitting solutions (wavelength channels)
    input_file_list = [f"{DATA_DIR}/{file}" for file in os.listdir(DATA_DIR)]

    # Initiate collection of white-noise (expected) and fit rms
    stderr_coll = []
    rms_coll = []
    bins = None

    # Loop over all files
    for file in input_file_list:
        residuals = lc_data(file)
        rms, bins, stderr = allan_values(residuals, binstep=10)

        # Add to collections
        stderr_coll.append(stderr)
        rms_coll.append(rms)

    # Plot normalized variances
    allan_plot(bins, rms_coll, stderr_coll)
    white_noise_plot(stderr_coll, bins)


def lc_data(filename):
    """Extract residuals from Eureka! LC fits"""
    data = pd.read_csv(filename, sep=" ", comment="#")

    return data["residuals"]


def allan_values(residuals, binstep=1, maxbins=None):
    """
    Compute the RMS, expected white-noise and bin intervals for a
    specified number of bin sizes and bin steps.
    """
    if maxbins is None:
        maxbins = residuals.shape[0] / 10

    binsz = np.arange(1, maxbins + binstep, step=binstep, dtype=int)
    nbins = np.zeros(binsz.size, dtype=int)
    rms = np.zeros(binsz.size)

    # TODO: This is straight-up stolen from the EUREKA! source code<
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


def allan_plot(bins, rms_list, stderr_list):
    """
    Generate combined Allan-Variance-Plot, where the x-axis represents
    the bin size in time and the y-axis the normalized RMS (to the
    expected white noise)
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot actual RMS, normalized to white-noise rms
    for idx in range(len(rms_list)):
        rms_norm = rms_list[idx] / stderr_list[idx][0]

        if idx == 0:
            ax.plot(bins, rms_norm, ls="-", lw=2, c="black", alpha=0.5,
                    label="RMS")
        else:
            ax.plot(bins, rms_norm, ls="-", lw=2, c="black", alpha=0.5)

    # Plot white noise expected RMS (identical because of normalization)
    ref_stderr = stderr_list[0]
    ax.plot(bins, ref_stderr / ref_stderr[0], ls="--", lw=2.5, c="tab:red",
            label="White-noise")

    ax.set(
        xscale="log", xlabel="Bin size (# of integrations)",
        yscale="log", ylabel="RMS"
    )
    plt.legend()
    plt.tight_layout()

    plt.savefig(f"{PLOT_DIR}/comb_allan_variance.{PLOT_TYPE}", dpi=600)


def white_noise_plot(stderr_coll, bins):
    """
    Plot the expected white-noise for all wavelength channels, color-coded
    by the channel number.
    """
    # Prepare color map by sampling with number of wavelength channels
    fig, ax = plt.subplots(figsize=(10, 8))
    color_map = mpl.colormaps['plasma'].resampled(len(stderr_coll))

    for idx in range(len(stderr_coll)):
        ax.plot(bins, stderr_coll[idx], lw=2,
                color=color_map.colors[idx])

    ax.set(
        xscale="log", xlabel="Bin size (# of integrations)",
        yscale="log", ylabel="Expected RMS (by channel)"
    )
    plt.tight_layout()

    # TODO: Colour bar indicating channel number
    plt.savefig(f"{PLOT_DIR}/comb_whitenoise.{PLOT_TYPE}", dpi=600)


if __name__ == "__main__":
    rc_setup()
    main()
