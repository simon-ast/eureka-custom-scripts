"""
OVERVIEW:

TBD

"""
import h5py as h5
import shutil
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from typing import Union
from util import rc_setup

# GLOBALS
PLOT_DIR = "plots/stage3_optspex"
FILE = "data/S3_wasp39b_ap4_bg8_FluxData_seg0000.h5"

# PARAMETERS
STD_TH = 10
SRC_POS = 13
SRC_HW = 4
SVAR_LIM = 21500


def main():
    # Extracts main necessary information from S3 data file
    stdspec, optspec, stdevs, optspec_rev, aux_spatial = \
        data_gen(SRC_POS, SRC_HW, FILE)

    # TODO: For now, I am trying to identify outliers by eye and mark
    #  them in the plot to correlate that
    bad_cols = np.genfromtxt("output/bad_indices_x.dat", delimiter=",",
                            dtype=int)

    # Prepare save directories for plots
    prep_save_dir(PLOT_DIR, bad_cols)

    for idx in bad_cols:
        # Plot the full standard and optimised LC for reference
        plot_lc_comparison(stdspec, optspec, idx)

        # Plot the spatial profile for the integration number range in
        # question (top-dictionary index by integration number!)
        for int_nr in aux_spatial.keys():
            spatial_nom = aux_spatial[int_nr]["spatial_nomask"]
            spatial_m = aux_spatial[int_nr]["spatial_masked"]
            plot_spatial_profile(spatial_m, spatial_nom, idx, int_nr)

        # Plot the stdev-variation with integration number
        stdev_col = stdevs[:, :, idx]
        int_range = np.arange(SVAR_LIM)
        plot_s_variation(stdev_col, int_range, idx, STD_TH)


def data_gen(src_row: int, src_hw: int, file_name: str):
    """
    Reproduces optimal spectral extraction from Eureka! and stores
    auxiliary data
    """
    # Reading in the data from S3, where "optspex" was already applied
    ap_lo, ap_hi = ap_range(src_row, src_hw)
    data = h5.File(file_name)

    # AUXILIARY DATA (CUTS ARE TO AP_HI + 1 BC OF LIST COMPREHENSION)
    # Electron-flux per pixel (TIME, Y, X)
    flux = data['flux'][()][:, ap_lo:ap_hi + 1, :]
    # Number of integrations
    n_int = flux.shape[0]
    # Background (TIME, Y, X)
    bg = data['bg'][()][:, ap_lo:ap_hi + 1, :]
    # Variance = noise^2 (TIME, Y, X)
    v0 = data['v0'][()][:, ap_lo:ap_hi + 1, :]
    # DQ mask from previous steps
    mask = data['mask'][()][:, ap_lo:ap_hi + 1, :]
    gain = 1

    # Median Frame only for source
    medflux = data['medflux'][()][ap_lo:ap_hi + 1, :]

    # PREVIOUSLY USED DATA FOR COMPARISON
    # Standard spectrum
    spectrum = data['stdspec'][()]
    # Optimal spectrum for comparison
    optspec = data['optspec'][()]

    # Calculate normalized median profile
    meddata = profile_meddata(medflux)

    # This is the calculation routine for optimal spectral extraction
    # used in Eureka!
    expected = np.array([meddata * spectrum[i] for i in range(n_int)])
    variance = np.ma.abs(expected + bg) / gain + v0
    stdevs = np.abs(flux - expected) * mask / np.sqrt(variance)

    # y-index of worst stdevs outlier in each column
    worst_loc = np.ma.argmax(stdevs, axis=1)

    # FOR EACH INTEGRATION
    # Determine a threshold, add to the mask and extract optimal spectrum
    # TODO: Right now, this is done 21500 times! by default
    n_int = 20

    # Initialise array for optimal spectrum and spatial data
    optspec_rev = []
    spatial_dict = {}

    # Loop over all time steps
    for int_nr in range(n_int):
        # Append the outlier mask
        outlier_mask(worst_loc[int_nr], stdevs[int_nr], mask[int_nr], STD_TH)

        # Use optimal spectral extraction step from Eureka!
        int_dict, optspex_step = optimal_spectral_extraction(
            meddata, int_nr, mask, flux, variance
        )

        # Append optimal spectral extraction array and spatial data
        optspec_rev.append([optspex_step])
        spatial_dict[f"{int_nr}"] = int_dict

    return spectrum, optspec, stdevs, optspec_rev, spatial_dict


def ap_range(source_row: int, source_hw: int):
    """Simple conversion from centre + half-width to absolute bound"""
    aperture_lo = source_row - source_hw
    aperture_hi = source_row + source_hw

    return aperture_lo, aperture_hi


def outlier_mask(pos_max_out: np.ndarray, sigma_array: np.ndarray,
                 masking_array: np.ndarray, threshold: int):
    """Masking outliers based on maximum sigma value"""
    # Loop over all "columns" (meaning x-positions)
    for id_x in range(pos_max_out.shape[0]):

        # Check of the standard deviation of the worst outlier is larger
        # than the threshold
        if sigma_array[pos_max_out[id_x], id_x] > threshold:
            # Add to the existing masking array
            masking_array[pos_max_out[id_x], id_x] = 0


def optimal_spectral_extraction(spatial_profile: np.ndarray, frame_nr: int,
                                masking_array: np.ndarray, flux: np.ndarray,
                                variance: np.ndarray):
    """DOC! REFERENCE TO EUREKA CODE!!"""
    # NO MASK
    temp_no = spatial_profile * flux[frame_nr] / variance[frame_nr]

    # WITH MASK
    temp = temp_no * masking_array[frame_nr]

    # DENOMINATOR
    temp_denom = spatial_profile * spatial_profile * masking_array[frame_nr] \
                 / variance[frame_nr]
    denom = np.sum(temp_denom, axis=0)

    # Optimal spectrum calculation (follows Eureka code)
    optimal_spectrum = np.sum(temp, axis=0) / denom

    # Prepare dictionary for later plots
    integration_dictionary = {
        "spatial_masked": temp,
        "spatial_nomask": temp_no
    }

    return integration_dictionary, optimal_spectrum


def profile_meddata(meddata):
    """
    Construct normalized spatial profile using median of all data frames.
    (SIMON): THIS IS TAKEN FROM "optspex.py" OF EUREKA FOR CONSISTENCY.
    """
    profile = np.copy(meddata)
    # Enforce positivity
    profile[np.where(profile < 0)] = 0
    # Normalize along spatial direction
    with np.errstate(divide='ignore', invalid='ignore'):
        profile /= np.sum(profile, axis=0)

    return profile


def prep_save_dir(root_dir, column_numbers):
    """Clean and repopulate save directory"""
    # List all existing folders
    previous_content = os.listdir(root_dir)

    # Remove all previous folders (this just does nothing if the
    # directory is empty)
    for folder in previous_content:
        shutil.rmtree(f"{root_dir}/{folder}")

    # Create new folders for each column in index list
    for column_idx in column_numbers:
        os.makedirs(f"{PLOT_DIR}/col{column_idx}")


def plot_s_variation(stdev_values: np.ndarray,
                     int_range: np.ndarray, col_idx: int,
                     sigma_threshold: Union[int, float]):
    """DOC!"""
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot the sigma-threshold used for masking
    ax.axhline(y=sigma_threshold, c="black", ls="--", lw=2)

    # Determine number of pixels in aperture
    ap_length = stdev_values.shape[1]

    # Generate custom colour map for each aperture pixel
    colours = mpl.colormaps['hsv'].resampled(ap_length)(range(ap_length))

    for aperture_px in range(ap_length):
        ax.plot(
            int_range, stdev_values[:, aperture_px][int_range],
            color=colours[aperture_px], label=f"{aperture_px}",
            marker="o", ls="-", lw=2
        )

    ax.set(
        xlabel="Integration number", ylabel="$S$-value (REF)",
        title=f"S-variation in aperture pixels of column {col_idx}",
        ylim=(0, sigma_threshold * 1.5)
    )

    mpl.rcParams["legend.frameon"] = "True"
    plt.legend(title="Aperture Pixel", ncols=ap_length // 3)
    mpl.rcParams["legend.frameon"] = "False"

    plt.tight_layout()

    plot_dir = f"{PLOT_DIR}/col{col_idx}"
    plt.savefig(f"{plot_dir}/col{col_idx}_s-variation.png", dpi=600)
    plt.close()


def plot_spatial_profile(masked, not_masked, bad_col, int_nr):
    """DOC"""
    rel_pxidx = np.arange(masked.shape[0])

    fig, ax = plt.subplots(figsize=(5, 4))

    ax.plot(rel_pxidx, not_masked[:, bad_col], color="tab:green",
            marker="o", label="Not masked")
    ax.plot(rel_pxidx, masked[:, bad_col], color="tab:red",
            marker="o", label="Masked")

    ax.set(
        xlabel="Aperture pixel index", ylabel="Arbitrary",
        title=f"Spatial profile of column {bad_col} at integration {int_nr}"
    )

    plt.legend()
    plt.tight_layout()

    plot_dir = f"{PLOT_DIR}/col{bad_col}"
    plt.savefig(f"{plot_dir}/col{bad_col}_int{int_nr}_spatial.png", dpi=600)
    plt.close()


def plot_lc_comparison(normal_spec: np.ndarray, opt_spec: np.ndarray,
                       col_idx: int):
    """Plot the full standard and optimised LC for a specified column"""
    # SANITY CHECK
    assert normal_spec.shape == opt_spec.shape, \
        "LC COMPARISON: Spectra do not have the same shape"

    # Standard parameters + flux selection and normalization
    integration_array = np.arange(normal_spec.shape[0])
    norm_stdspec = normal_spec[:, col_idx] / np.mean(normal_spec[:, col_idx])
    norm_optspec = opt_spec[:, col_idx] / np.mean(opt_spec[:, col_idx])

    # PLOTTING ROUTINE
    fig, ax = plt.subplots(figsize=(8, 5))

    # Plot standard spectrum
    ax.scatter(integration_array, norm_stdspec, marker="o", color="tab:blue",
               alpha=0.6, zorder=0, label="Standard spec.")

    # Plot optimal spectral extraction from Eureka!
    ax.scatter(integration_array, norm_optspec, marker="^", color="tab:orange",
               alpha=0.6, zorder=1, label="Optimal spec.")

    # Auxiliary labels
    ax.set(
        xlabel="Integration number", ylabel="Normalised Flux",
        title=f"2D LC for column {col_idx}"
    )

    plt.legend()
    plt.tight_layout()

    plot_dir = f"{PLOT_DIR}/col{col_idx}"
    plt.savefig(f"{plot_dir}/col{col_idx}_2d-LC.png", dpi=600)
    plt.close()


if __name__ == "__main__":
    rc_setup()
    main()

