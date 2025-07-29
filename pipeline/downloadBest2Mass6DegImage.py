import os
import configparser
from astropy.wcs import WCS
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.units as u
from astropy.utils.data import download_file
import numpy as np
from datetime import datetime, timezone, UTC
from dateutil import tz
import time
import pyvo as vo
import re


to_zone = tz.gettz('America/Los_Angeles')

import modules.utils.gdip_subs as util
import pipeline.gdip_referenceImageSubs as rfis
import pipeline.gdip_differenceImageSubs as dfis

start_time_benchmark = time.time()
start_time_benchmark_at_start = start_time_benchmark


swname = "generalDifferenceImagePipeline.py"
swvers = "1.0"
cfg_filename_only = "generalDifferenceImagePipeline.ini"

print("swname =", swname)
print("swvers =", swvers)
print("cfg_filename_only =", cfg_filename_only)


# Assume the RAPID-pipeline git repo is located under /code.

rapid_sw = "/code"


# Specify science image.

filename_science_image = 'ADP.2022-07-27T14_56_30.297.fits'
hdu_index_science = 1                                                    # Second HDU
jid = 1


# Compute processing datetime (UT) and processing datetime (Pacific time).

datetime_utc_now = datetime.now(UTC)
proc_utc_datetime = datetime_utc_now.strftime('%Y-%m-%dT%H:%M:%SZ')
datetime_pt_now = datetime_utc_now.replace(tzinfo=timezone.utc).astimezone(tz=to_zone)
proc_pt_datetime_started = datetime_pt_now.strftime('%Y-%m-%dT%H:%M:%S PT')

print("proc_utc_datetime =",proc_utc_datetime)
print("proc_pt_datetime_started =",proc_pt_datetime_started)


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


    # Read FITS file

    with fits.open(filename_science_image) as hdul:

        filter_science = hdul[0].header["FILTER"].strip()

        print("filter_science =",filter_science)

        hdr = hdul[hdu_index_science].header

        radecsys = hdr['RADECSYS']
        print("radecsys =",radecsys)
        hdr.remove('RADECSYS', remove_all=True)
        hdr['RADESYSA'] = radecsys

        hdr.remove('PROJP1', remove_all=True)
        hdr.remove('PROJP3', remove_all=True)
        hdr.remove('PROJP5', remove_all=True)

        w_sci = WCS(hdr) # Initialize WCS object from FITS header

    print(w_sci)

    print("CTYPE = ",w_sci.wcs.crpix)

    naxis1 = hdr['NAXIS1']
    naxis2 = hdr['NAXIS2']

    print("naxis1,naxis2 =",naxis1,naxis2)

    crpix1 = w_sci.wcs.crpix[0]
    crpix2 = w_sci.wcs.crpix[1]


    # Example of converting pixel coordinates to celestial coordinates
    # The following should reproduce CRVAL1,CRVAL2.

    pixel_x, pixel_y = crpix1 - 1, crpix2 - 1
    celestial_coords = w_sci.pixel_to_world(pixel_x, pixel_y)
    print(f"CRVAL1,CRVAL2 Pixel ({pixel_x}, {pixel_y}) corresponds to {celestial_coords.ra.deg:.12f} RA and {celestial_coords.dec.deg:.12f} Dec.")


    # Compute pixel coordinates of science-image center and four corners.

    x0,y0,x1,y1,x2,y2,x3,y3,x4,y4 = util.compute_pix_image_center_and_four_corners(naxis1,naxis2)


    # Compute RA,Dec of center and four corners of science image.

    (ra0,dec0,ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4) = util.compute_sky_image_center_and_four_corners(w_sci,x0,y0,x1,y1,x2,y2,x3,y3,x4,y4)


    # Specify sky position of interest and desired filter.

    ra = ra0
    dec = dec0
    pos = SkyCoord(ra=ra, dec=dec, unit='deg')
    print("RA,Dec of image center =",ra,dec)

    # Query for 2MASS mosaic_sixdeg images that overlap sky position within 1.0 arcseconds.

    download_url = "https://irsa.ipac.caltech.edu/SIA?COLLECTION=twomass_sixdeg&POS=circle+85.53725422808287+-68.72445403392891+0.0003&RESPONSEFORMAT=IPAC_TABLE"

    table_fits_file_ref = "twomass_sixdeg.ipac"
    curl_cmd = "curl --output " + table_fits_file_ref + " \"" + download_url + "\""
    print("curl_cmd =",curl_cmd)

    return_code = os.system(curl_cmd)
    print(f"Command exited with code: {return_code}")

    im_table = ascii.read(table_fits_file_ref, format='ipac')

    print(im_table)


    # Columns of interest in table: access_url and energy_bandpassname.

    # Get the returned 2MASS image(s) for the filter of interest.
    # Pick the one that overlaps the science image the most.

    max_percent_overlap_area = 0.0
    intersecting_polygon_file = "intersecting_polygon.txt"

    for i in range(len(im_table)):

        if "fits" not in im_table[i]['access_url']:
            continue

        if im_table[i]['energy_bandpassname'].lower() != filter_science.lower():
            continue

        download_url = im_table[i]['access_url']
        print("download_url =",download_url)

        filename_match = re.match(r".+/(.+?.fits)", download_url)

        try:
            provisional_fits_file_ref = filename_match.group(1)
            print("provisional_fits_file_ref =",provisional_fits_file_ref)

            curl_cmd = "curl --output " + provisional_fits_file_ref + " \"" + download_url + "\""
            print("curl_cmd =",curl_cmd)

            return_code = os.system(curl_cmd)
            print(f"Command exited with code: {return_code}")

        except:
            print("*** Error: No FITS filename match found; quitting...")
            exit(64)


        # Take the absolute value in case the reference image is flipped and polygon vertex order is reversed.

        percent_overlap_area = util.compute_image_overlap_area(w_sci,naxis1,naxis2,x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,provisional_fits_file_ref)

        if percent_overlap_area > max_percent_overlap_area:

            try:
                delete_cmd = "rm -f " + fits_file_ref
                print("delete_cmd =",delete_cmd)

                return_code = os.system(delete_cmd)
                print(f"Command exited with code: {return_code}")
            except:
                pass

            max_percent_overlap_area = percent_overlap_area
            fits_file_ref = provisional_fits_file_ref

            copy_cmd = "cp -f outfile intersecting_polygon_file"
            print("copy_cmd =",copy_cmd)

            return_code = os.system(copy_cmd)
            print(f"Command exited with code: {return_code}")


        else:
            delete_cmd = "rm -f " + provisional_fits_file_ref
            print("delete_cmd =",delete_cmd)

            return_code = os.system(delete_cmd)
            print(f"Command exited with code: {return_code}")


    # The reference image with maximal overlap area and matching filter/band has been selected.

    print("\n")
    print("Reference image =",fits_file_ref)
    print(f"max_percent_overlap_area = {max_percent_overlap_area:.2f}")

