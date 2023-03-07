"""
OVERVIEW:

Generates a precision plot showing MAD values for each pixel column

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from util import rc_setup

# GLOBALS
PLOT_DIR = "plots/"
PLOT_TYPE = "png"
MAD_FILE = "output/mad_values.dat"
ANNOTATE = False


def main():
    """DOC!"""
    # Generate data frame with wavelength and MAD values
    mad_data = read_mad(MAD_FILE)

    # Calculate precision
    fit_par, _ = fit_precision(mad_data.wavel, mad_data.mad_values)

    # Identify outliers from the fitted precision
    # TODO: For now, I am trying to identify outliers by eye and mark
    #  them in the plot to correlate that
    bad_idx = np.genfromtxt("output/bad_indices_x.dat", delimiter=",",
                            dtype=int)

    # Plot a tidy plot
    plot_precision(mad_data.wavel, mad_data.mad_values, fit_par, bad_idx)


def read_mad(filename):
    """Read in MAD data file"""
    # Read in the data file
    mad_data = pd.read_csv(filename, sep="\t")

    # Rename the columns for easier reference
    mad_data.columns = ["wavel", "mad_values"]

    return mad_data


def fit_precision(wavelength, mad_values):
    """Not in use for now"""
    # TODO: Find an automated way of flagging and cross-checking (!) bad pixels

    # Fit an exponential with scipy
    res = curve_fit(exponential, wavelength, mad_values)

    popt = res[0]
    pcov = res[1]

    return popt, pcov


def large_outlier_rejection(mad_values, threshold):
    """Not in use for NOW"""
    n_dp = len(mad_values)

    temp_1 = mad_values
    temp_2 = np.roll(mad_values, -1)

    diff_up = np.abs((temp_1 - temp_2))[:-1]
    diff_do = np.abs((temp_2 - temp_1))[1:]
    #print(len(diff_up), len(diff_do))

    test_1 = np.where(diff_up > threshold)
    test_2 = np.where(diff_do > threshold)

    test_3 = np.in1d(test_1, test_2)
    eval = test_3 * test_1[0]
    test_4 = eval[eval != 0] - 1

    #plt.scatter(np.arange(n_dp - 1), diff_up)
    #plt.scatter(np.arange(n_dp - 1), diff_do)
    #plt.show()

    return test_4


def plot_precision(wavelength, mad_values, fit_par, flagged_indices):
    """
    Plots the MAD value for each light curve, and marks pixels that
    represent large outliers
    """
    fig, ax = plt.subplots(figsize=(6, 3))

    # Plotting parameters
    mad_ms = 15

    # Plot MAD values
    ax.scatter(wavelength, mad_values / 1e3, color="tab:blue", s=mad_ms)
    ax.set(xlabel="Wavelength [$\\mu \\mathrm{m}$]",
           ylabel="MAD [ppt]", ylim=(2, 55))

    # Overplot flagged indices
    ax.scatter(wavelength[flagged_indices], mad_values[flagged_indices] / 1e3,
               c="tab:red", s=mad_ms, marker="x", label="Flagged Values")

    # TODO: Use annotations if necessary to find bad values
    if ANNOTATE is True:
        for i in range(len(wavelength)):
            ax.annotate(text=f"{i}", xy=(wavelength[i], mad_values[i] / 1e3))

    # Plot fitted precision (TODO: Excluded for now)
    # wl_prec = exponential(wavelength, *fit_par) / 1e3
    # ax.plot(wavelength, wl_prec, ls="--", lw=2.5, c="tab:red")

    plt.legend()
    plt.tight_layout()
    plot_loc = f"{PLOT_DIR}/stage4_precision"
    plt.savefig(f"{plot_loc}/eucus_precision.{PLOT_TYPE}", dpi=600)


def exponential(x, a, b, c):
    """Function used for fitting (YANI; CITATION?)"""
    return a * np.exp(b * x)


if __name__ == "__main__":
    # Custom plot adjustments from util.py
    rc_setup()

    # Call of main routine
    main()
