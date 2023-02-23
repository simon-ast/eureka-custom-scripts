"""
OVERVIEW:

Customised plotting routine for 2D (dynamic) light curves from the Stage 3
data products produced by 'Eureka!', and the possibility to plot individual
light curves. This routine reads from the
'S3*_SpecData.hdf5' files produced in that stage.

Also includes precision plot for MAD [ppm].

"""

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
from util import rc_setup

# DIRECTORIES RELATIVE TO "eureka_custom_plots"!
DATA_DIR = "data"
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
    wavel = data["wave_1d"][()]

    # Plot and save dynamic light curve and precision plot
    plot_dynamic_lc(wavel, n_int, normspec, "RdYlBu_r", "wavelength")
    plot_precision(wavel, normspec)

    # Plot 1D light curves for specified pixel column
    idx_bad = [7, 12]
    for idx in idx_bad:
        plot_1d_lightcurve(n_int, normspec, idx)


def plot_dynamic_lc(x_array, y_array, z_array,
                    colormap, save_append, x_limits=None):
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


def plot_precision(wavelength, norm_flux):
    """Plots MAD for each time series column in the 2d LC."""
    fig, ax = plt.subplots(figsize=(6, 3))

    # Basic parameters
    n_int, n_pix = norm_flux.shape

    # Find wavelength cutoff at 5.3 microns (don't forget that list
    # slicing is exclusive for the stopping index!)
    end_idx = np.where(wavelength > 5.3)[0][0]
    wavelength = wavelength[:end_idx]

    # Deviation from next neighbour of every point in every column
    # Restricted by wavelength cutoff
    ediff1d = np.array([
        np.ediff1d(norm_flux[:, i]) for i in range(n_pix)
    ])[:end_idx]

    # Median Absolute Deviation in PPM
    mad_array = np.ma.median(np.abs(ediff1d), axis=1) * 1e6

    # Plot everything nice and tidy
    ax.scatter(wavelength, mad_array / 1e3, color="tab:blue", s=15)
    ax.set(xlabel="Wavelength [$\\mu \\mathrm{m}$]",
           ylabel="MAD [ppt]")

    plt.tight_layout()
    plot_loc = f"{PLOT_DIR}/stage3_2d-lightcurve"
    plt.savefig(f"{plot_loc}/eucus_precision.{PLOT_TYPE}", dpi=600)


def plot_1d_lightcurve(n_int, spec, px_idx):
    """Plot 1D data set from the 2D LC"""
    fig, ax = plt.subplots(figsize=(10, 3))
    light_c = spec[:, px_idx]

    ax.scatter(n_int, light_c, alpha=0.5)
    ax.set(
        title=f"Relative pixel number {px_idx + 1} (Dispersion direction)",
        xlabel="Integration number",
        ylabel="Normalized Flux", ylim=(0.8, 1.2)
    )

    plt.tight_layout()
    plot_loc = f"{PLOT_DIR}/stage3_2d-lightcurve"
    plt.savefig(f"{plot_loc}/eucus_1dlc_px{px_idx}.{PLOT_TYPE}", dpi=600)


if __name__ == "__main__":
    # Custom plot adjustments from util.py
    rc_setup()

    # Call of main routine
    main()
