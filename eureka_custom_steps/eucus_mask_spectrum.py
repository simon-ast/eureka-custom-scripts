from astropy.io import fits
import numpy as np

# GLOBALS
FILE_DIR = "output"


def alter_fits_file(filename, mask_array):
    """DOCSTRING!"""
    with fits.open(filename, mode='update') as hdul:
        data = hdul[1].data
        for nint in range(data.shape[0]):
            for ngroup in range(data.shape[1]):
                mask_2d_spectra(data[nint][ngroup], mask_array)
        hdul.flush()


def mask_2d_spectra(data, mask_array):
    """DOC!"""
    for entry in mask_array:
        # Subtract 1 from mask entries! Pixel coordinates start at 1,
        # array entries start at 0
        x, y = entry[0] - 1, entry[1] - 1
        data[y][x] = 0


def main():
    file_list = [
        "jw01366004001_04101_00001-seg001_nrs1_uncal",
        "jw01366004001_04101_00001-seg002_nrs1_uncal",
        "jw01366004001_04101_00001-seg003_nrs1_uncal",
        "jw01366004001_04101_00001-seg004_nrs1_uncal"
    ]
    mask = [
        [57, 27], [72, 4], [106, 13], [121, 19], [131, 22], [159, 29],
        [196, 3], [201, 28], [290, 25], [340, 7], [379, 11], [448, 32],
        [447, 27], [478, 25], [502, 11], [505, 22]
    ]

    for file in file_list:
        file_name = f"{file}.fits"
        alter_fits_file(file_name, mask)
    pass


if __name__ == "__main__":
    main()
