"""
OVERVIEW:

Customised plotting routine for 2D (dynamic) light curves from the Stage 3
data products produced by 'Eureka!'. This routine reads from the
'S3*_SpecData.hdf5' files produced in that stage.

Also calculates and saves MAD values for each light curve.

"""

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from util import rc_setup

# DIRECTORIES RELATIVE TO "eureka_custom_plots"!
DATA_DIR = "data"
OUT_DIR = "output"
FILE_NAME = "S3_wasp39b_ap4_bg8_SpecData.h5"
PLOT_DIR = "plots/"
PLOT_TYPE = "png"


def main():
    # Complete data frame
    data = h5.File(f"{DATA_DIR}/{FILE_NAME}")

    # SPECTRUM AFTER STEP "optimal spectral extraction"
    # Has the shape (TIME, DISPERSION-DIRECTION)
    spec = data["optspec"][()]

    # Normalization for each time-column
    normspec = spec / np.ma.mean(spec, axis=0)

    # AUXILIARY DATA: Number of integrations, x-direction pixel-number,
    # assigned wavelength in microns
    n_int = np.arange(spec.shape[0])
    n_pix = np.arange(spec.shape[1]) + 1    # Should start at 1
    wavel = data["wave_1d"][()]

    # TODO: Read-in bad column indices
    bad_col = np.genfromtxt("output/bad_indices_x.dat", delimiter=",",
                            dtype=int)

    # Plot and save dynamic light curve and precision plot
    plot_dynamic_lc(wavel, n_int, normspec, "RdYlBu_r", "wavelength", bad_col)

    # Calculate MAD values for each normalised light curve
    save_mad(wavel, None, 5.3, normspec, "mad_values")


def plot_dynamic_lc(x_array, y_array, z_array,
                    colormap, save_append,
                    bad_px_idx, x_limits=None):
    """
    Plotting routine for dynamic light curve (could be either with
    respect to wavelength or pixel position.
    """
    # PLOTTING 2D LIGHT CURVE, can be cut into custom segments
    fig, ax = plt.subplots(figsize=(12, 6))
    dyn_lc = ax.pcolormesh(
        x_array, y_array, z_array,
        cmap=colormap, vmin=0.97, vmax=1.03
    )

    # TODO: Plot location (x) of flagged pixel columns
    for idx in bad_px_idx:
        ax.axvline(x_array[idx], c="black", alpha=0.5, zorder=10000)

    plt.colorbar(
        dyn_lc, label="Flux normalized to mean per channel"
    )

    # Set x-limits corresponding to either wavelength or pixel-number
    ax.set(xlim=x_limits, ylabel="Integration number")
    if save_append == "wavelength":
        ax.set(xlabel="Wavelength $[\\mu \\mathrm{m}]$")
    elif save_append == "pixel":
        ax.set(xlabel="Pixel position of source")

    plt.tight_layout()

    # Plots are saved in "plots/stage3_2d-lightcurve"
    plot_loc = f"{PLOT_DIR}/stage3_2d-lightcurve"
    plt.savefig(
        f"{plot_loc}/eucus_2Dlightcurve_{save_append}.{PLOT_TYPE}",
        dpi=600
    )


def save_mad(wavelength, wl_lo, wl_hi, normalised_flux, output_name):
    """
    Calculates the Median Absolute Deviation for each spectral light 
    curve and saves the results to a file
    """
    # Basic parameters:
    # number of integrations and pixel dimension in dispersion direction
    n_int, n_pix = normalised_flux.shape

    # Find wavelength cutoff at 5.3 microns (don't forget that list
    # slicing is exclusive for the stopping index!)
    start_idx, end_idx = wavelength_cutoffs(wavelength, wl_lo, wl_hi)
    wavelength = wavelength[:end_idx]

    # Deviation from next neighbour of every point in every column
    # Restricted by wavelength cutoff
    ediff1d = np.array([
        np.ediff1d(normalised_flux[:, i]) for i in range(n_pix)
    ])[:end_idx]

    # Median Absolute Deviation in PPM
    mad_array = np.ma.median(np.abs(ediff1d), axis=1) * 1e6

    # Wavelengths and corresponding MAD values to file
    cols = ["wl [micr]", "MAD [ppm]"]
    stacked_data = np.stack((wavelength, mad_array), axis=-1)
    data_frame = pd.DataFrame(stacked_data, columns=cols)
    data_frame.to_csv(f"{OUT_DIR}/{output_name}.dat", sep="\t", index=False)


def wavelength_cutoffs(wavelength, lower_bound, upper_bound):
    """Calculates array indices for upper and lower wavelength boundaries"""
    wavebound_low = lower_bound
    if lower_bound is not None:
        wavebound_low = np.where(wavelength < lower_bound)[0][0]

    wavebound_high = upper_bound
    if upper_bound is not None:
        wavebound_high = np.where(wavelength > upper_bound)[0][0]

    return wavebound_low, wavebound_high


if __name__ == "__main__":
    # Custom plot adjustments from util.py
    rc_setup()

    # Call of main routine
    main()
