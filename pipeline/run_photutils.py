import os
import configparser
from astropy.wcs import WCS
from astropy.io import fits
from astropy.io import ascii
from datetime import datetime, timezone, UTC
from dateutil import tz
import time


to_zone = tz.gettz('America/Los_Angeles')

import modules.utils.gdip_subs as util


start_time_benchmark = time.time()
start_time_benchmark_at_start = start_time_benchmark


swname = "run_photutils.py"
swvers = "1.0"
cfg_filename_only = "generalDifferenceImagePipeline.ini"

print("swname =", swname)
print("swvers =", swvers)
print("cfg_filename_only =", cfg_filename_only)


# Get the configuration parameters.

diffimg_sw = os.getenv('DIFFIMG_SW')

if diffimg_sw is None:

    print("*** Error: Env. var. DIFFIMG_SW not set; quitting...")
    exit(64)

diffimg_work = os.getenv('DIFFIMG_WORK')

if diffimg_work is None:

    print("*** Error: Env. var. DIFFIMG_WORK not set; quitting...")
    exit(64)

cfg_path = diffimg_sw + "/cdf"

print("diffimg_sw =",diffimg_sw)
print("cfg_path =",cfg_path)


# Read input parameters from .ini file.

config_input_filename = cfg_path + "/" + cfg_filename_only
config_input = configparser.ConfigParser()
config_input.read(config_input_filename)

print("config_input_filename =",config_input_filename)


#-------------------------------------------------------------------------------------------------------------
# Main program.
#-------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':


    # Generate PSF-fit catalog for ZOGY difference image using photutils.  No background subtraction is done.

    filename_bkg_subbed_science_image = config_input['BKGEST']['filename_bkg_subbed_science_image']

    psfcat_diffimage_dict = config_input['PSFCAT_DIFFIMAGE']

    n_clip_sigma = float(psfcat_diffimage_dict["n_clip_sigma"])
    n_thresh_sigma = float(psfcat_diffimage_dict["n_thresh_sigma"])

    fwhm = float(psfcat_diffimage_dict["fwhm"])
    fit_shape_str = psfcat_diffimage_dict["fit_shape"]
    fit_shape = tuple(int(x) for x in fit_shape_str.replace("(","").replace(")","").replace(" ", "").split(','))
    aperture_radius = float(psfcat_diffimage_dict["aperture_radius"])


    input_img_filename = psfcat_diffimage_dict["input_img_filename"]
    input_unc_filename = psfcat_diffimage_dict["input_unc_filename"]
    input_psf_filename = psfcat_diffimage_dict["input_psf_filename"]
    output_psfcat_filename = psfcat_diffimage_dict["output_psfcat_filename"]
    output_psfcat_finder_filename = psfcat_diffimage_dict["output_psfcat_finder_filename"]
    output_psfcat_residual_filename = psfcat_diffimage_dict["output_psfcat_residual_filename"]

    psfcat_flag,phot,psfphot = util.compute_diffimage_psf_catalog(n_clip_sigma,
                                                                  n_thresh_sigma,
                                                                  fwhm,
                                                                  fit_shape,
                                                                  aperture_radius,
                                                                  input_img_filename,
                                                                  input_unc_filename,
                                                                  input_psf_filename,
                                                                  output_psfcat_residual_filename)

    print("psfcat_flag =",psfcat_flag)

    if psfcat_flag:


        # Output psf-fit catalog is an PSFPhotometry astropy table with the PSF-fitting results
        # merged with the DAOStarFinder astropy table.
        # Output columns are documentated at
        # https://photutils.readthedocs.io/en/latest/api/photutils.psf.PSFPhotometry.html
        # https://photutils.readthedocs.io/en/stable/api/photutils.detection.DAOStarFinder.html

        try:
            phot['x_init'].info.format = '.4f'
            phot['y_init'].info.format = '.4f'
            phot['flux_init'].info.format = '.6f'
            phot['flux_fit'].info.format = '.6f'
            phot['x_err'].info.format = '.4f'
            phot['y_err'].info.format = '.4f'
            phot['flux_err'].info.format = '.5f'
            phot['qfit'].info.format = '.4f'
            phot['cfit'].info.format = '.4f'

            print(phot[('id', 'x_fit', 'y_fit', 'flux_fit','x_err', 'y_err', 'flux_err', 'npixfit', 'qfit', 'cfit', 'flags')])


            # Compute sky coordinates for given pixel coordinates.

            ra,dec = util.computeSkyCoordsFromPixelCoords(filename_bkg_subbed_science_image,
                                                          list(phot['x_fit']),
                                                          list(phot['y_fit']))

            phot['x_fit'].info.format = '.4f'
            phot['y_fit'].info.format = '.4f'
            phot.add_column(ra, name='ra')
            phot.add_column(dec, name='dec')
            phot['ra'].info.format = '.6f'
            phot['dec'].info.format = '.6f'


            # Write PSF-fit photometry catalog in astropy table to text file.

            print("output_psfcat_filename = ", output_psfcat_filename)

            ascii.write(phot, output_psfcat_filename, overwrite=True)


            # Write PSF-fit finder catalog in astropy table to text file.

            print("output_psfcat_finder_filename = ", output_psfcat_finder_filename)

            ascii.write(psfphot.finder_results, output_psfcat_finder_filename, overwrite=True)

        except Exception as e:
            print(f"PSF-fit PSFPhotometry and DAOStarFinder catalogs: An unexpected error occurred: {e}")


    # Code-timing benchmark.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds after generating PSF-fit catalog on difference image =",
        end_time_benchmark - start_time_benchmark)
    start_time_benchmark = end_time_benchmark



    # Termination.

    terminating_exitcode = 0
    exit(terminating_exitcode)
