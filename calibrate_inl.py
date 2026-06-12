import numpy as np
import os
from astropy.io import fits
import argparse

def eval_fourier_counts(data, theta, periods):
    fs = np.zeros_like(data, dtype=float)
    K = periods.size
    
    for k in range(K):
        w = 2.0 * np.pi / periods[k]
        s_coef = theta[2*k + 0]
        c_coef = theta[2*k + 1]
        arg = w * data
        fs += s_coef * np.sin(arg) + c_coef * np.cos(arg)

    return fs


def main(args):
    input_path, output_path = args
    if not os.path.exists(output_path):
        os.makedirs(output_path)
            
    for fits_file in os.listdir(input_path):
        if fits_file.endswith('.fits'):
            data = fits.open(os.path.join(input_path, fits_file))
            uncal_data = data[1].data
            periods_grid = np.asarray([1024.0/3.0, 1024.0/2.0, 1024.0])
            theta = np.load("fourier_series_amplitudes.npy")
            correction = 1.0 + eval_fourier_counts(uncal_data, theta, periods_grid)
            data_corrected = uncal_data / correction
            
            # Save the corrected data to a new FITS file
            hdu = fits.PrimaryHDU(data_corrected)
            hdu.writeto(os.path.join(output_path, f'{fits_file}'), overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calibrate INL for FITS files.')
    parser.add_argument('input_path', type=str, help='Path to the input directory containing FITS files.')
    parser.add_argument('output_path', type=str, help='Path to the output directory for corrected FITS files.')
    
    args = parser.parse_args()
    main((args.input_path, args.output_path))
