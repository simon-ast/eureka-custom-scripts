"""
OVERVIEW:

- Perform weighted-average binning of transmission spectrum after LC
fitting has been applied.
- Plot transmission spectra comparison if desired

"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from util import rc_setup

# GLOBALS
DATA_DIR = "data/S6_spectral_binning"
PLOT_DIR = "plots/S6_spectral_binning"
OUT_DIR = "output/S6_spectral_binning"

# PARAMETERS
BIN_SIZE = 4
# TODO: READ BAD PIXEL COLUMNS FROM FILE
MASK_IDX = [7, 12, 14, 26, 40, 48, 69]


def main():
    # Determine files to be analysed
    optspex_file = f"{DATA_DIR}/eureka_optspex_nativeres.txt"
    stdspec_file = f"{DATA_DIR}/eureka_stdspec_nativeres.txt"

    # Read in transmission spectrum
    optspex_nat = read_spectrum(optspex_file)
    stdspec_nat = read_spectrum(stdspec_file)

    # TODO: Need to mask spectrum
    optspex_nat_ma = mask_data(optspex_nat, MASK_IDX)

    # Compute weighted arithmetic mean of spectra
    optspex_bin = spectrum_binning(optspex_nat_ma, bin_size=BIN_SIZE)
    stdspec_bin = spectrum_binning(stdspec_nat, bin_size=BIN_SIZE)

    # Plot the desired spectra and comparisons
    spectra = [optspex_nat, optspex_bin]
    labels = ["OS_NAT", "OS_BINT"]
    plot_spectrum_total(spectra, labels)
    plt.show()

    spectra = [stdspec_nat, stdspec_bin]
    labels = ["ST_NAT", "ST_BINT"]
    plot_spectrum_total(spectra, labels)
    plt.show()


def read_spectrum(filename):
    """Read data files that are results from Eureka!"""
    data_frame = pd.read_csv(filename, comment="#", delimiter=" ")

    return data_frame


def mask_data(spectrum, mask_array):
    """Mask specified indices in array as NaNs"""

    for ma_idx in mask_array:
        temp_data = spectrum.iloc[ma_idx]

        for col in temp_data.keys():
            temp_data.loc[col] = np.nan

    return spectrum


def spectrum_binning(pd_spectrum, bin_size):
    """Bin spectrum into determined blocks and calculate WAM"""
    # TODO: For now take average of errors, think of a better weighting
    #  mechanism!
    n_datapoints = pd_spectrum.shape[0]

    # Set-up empty pandas DF for binned data
    binned_spectrum = pd.DataFrame(columns=pd_spectrum.columns)

    idx = 0
    while idx <= n_datapoints - bin_size:
        # TODO: This might leave out some data points at the end
        temp_frame_raw = pd_spectrum[idx:idx + bin_size]
        temp_frame = temp_frame_raw.dropna().reset_index()
        temp_size = temp_frame.shape[0]

        # Calculate wavelength centre and bounds
        temp_wavel = np.sum(temp_frame.wavelength) / temp_size

        temp_wavel_lo = temp_frame.wavelength[0]
        temp_wavel_hi = temp_frame.wavelength[temp_size - 1]
        temp_binwidth = temp_wavel_hi - temp_wavel_lo

        # Calculate statistical parameters (weighted mean and std)
        temp_wam, temp_wam_std = weighted_arithmetic_mean(temp_frame)

        # Collect temporary data and construct temporary data frame
        temp_data = [temp_wavel, temp_binwidth,
                     temp_wam, temp_wam_std, temp_wam_std]
        temp_df = pd.DataFrame(data=[temp_data], columns=pd_spectrum.columns)

        # Append to total data frame
        binned_spectrum = pd.concat([binned_spectrum, temp_df])

        # Clean-up
        idx += bin_size

    return binned_spectrum


def weighted_arithmetic_mean(subframe):
    """Calculate WAM and its standard deviation"""
    # Store data in temporary variable
    data = subframe["rp^2_value"]

    # TODO: Temporary - variance from positive and negative error
    temp_var = ((subframe["rp^2_errorpos"] + subframe["rp^2_errorneg"])
                / 2) ** 2

    # Creating weights by inverse variance
    weights = 1 / temp_var

    # STD for WAM
    wam_std = np.sqrt(1 / np.sum(weights, axis=0))

    # Weighted arithmetic mean of subsets
    wam = np.sum(data * weights, axis=0) / np.sum(weights, axis=0)

    return wam, wam_std


def plot_spectrum_ind(axis, spectrum, label):
    """Plot individual spectrum onto axis object"""
    axis.errorbar(
        x=spectrum.wavelength, y=spectrum["rp^2_value"],
        xerr=[spectrum.bin_width, spectrum.bin_width],
        yerr=[spectrum["rp^2_errorneg"], spectrum["rp^2_errorpos"]],
        label=label, ls=" ", marker="o"
    )


def plot_spectrum_total(spectrum_df_arr, label_arr):
    """Plot all desired spectra in single plot"""
    # Number of spectra submitted
    n_spec = len(spectrum_df_arr)

    # Plot all spectra
    fig, ax = plt.subplots(figsize=(10, 6))

    for idx in range(n_spec):
        plot_spectrum_ind(ax, spectrum_df_arr[idx], label_arr[idx])

    ax.set(
        xlabel="$\\lambda$ [$\\mu \\mathrm{m}$]",
        ylabel="Transit Depth"
    )

    plt.legend()
    plt.tight_layout()


if __name__ == "__main__":
    rc_setup()
    main()
