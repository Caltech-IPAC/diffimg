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
from datetime import datetime, timezone
from dateutil import tz
import time
import pyvo as vo
import re


to_zone = tz.gettz('America/Los_Angeles')

import modules.utils.pipeline_subs as util
#import pipeline.referenceImageSubs as rfis
#import pipeline.differenceImageSubs as dfis

start_time_benchmark = time.time()
start_time_benchmark_at_start = start_time_benchmark


swname = "generalDifferenceImagePipeline.py"
swvers = "1.0"

print("swname =", swname)
print("swvers =", swvers)


# Compute processing datetime (UT) and processing datetime (Pacific time).

datetime_utc_now = datetime.utcnow()
proc_utc_datetime = datetime_utc_now.strftime('%Y-%m-%dT%H:%M:%SZ')
datetime_pt_now = datetime_utc_now.replace(tzinfo=timezone.utc).astimezone(tz=to_zone)
proc_pt_datetime_started = datetime_pt_now.strftime('%Y-%m-%dT%H:%M:%S PT')

print("proc_utc_datetime =",proc_utc_datetime)
print("proc_pt_datetime_started =",proc_pt_datetime_started)


# Specify science image.

input_fits_file = 'ADP.2022-07-27T14_56_30.297.fits'
hdr_num = 1                                                    # Second HDU




#-------------------------------------------------------------------------------------------------------------
# Main program.
#-------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':


    # Read FITS file

    with fits.open(input_fits_file) as hdul:

        filter_science = hdul[0].header["FILTER"].strip()

        print("filter_science =",filter_science)

        hdr = hdul[hdr_num].header

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


    # Query for 2MASS images that overlap sky position within 1.0 arcseconds.

    twomass_service = vo.dal.SIAService("https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?type=at&ds=asky&")
    im_table = twomass_service.search(pos=pos, size=1.0*u.arcsec)
    print(im_table.to_table())


    # Columns in table.
    # name,download,center_ra,center_dec,naxes,naxis,scale,format,crpix,crval,crota2,band,bref,bhi,blo,pers_art,glint_art,type,dataset,pixflags,id,scntr,date,hem,scan,image,ut_date,coadd_key,seesh,magzp,msnr10,bin


    # Get the first returned 2MASS image for the filter of interest.

    for i in range(len(im_table)):
        if im_table[i]['band'] == filter_science:
            break
    print(im_table[i].getdataurl())


    # Download the 2MASS image and decompress it.
    #
    # curl --output hi0600256.fits.gz "https://irsa.ipac.caltech.edu:443/cgi-bin/2MASS/IM/nph-im?ds=asky&atdir=/ti08&dh=000616n&scan=060&name=hi0600256.fits"
    # gunzip hi0600256.fits.gz

    download_url = im_table[i].getdataurl()

    filename_match = re.match(r".+?name\=(.+?.fits)", download_url)

    try:
        fits_file_ref = filename_match.group(1)
        print("fits_file_ref =",fits_file_ref)

        gz_fits_file_ref = fits_file_ref + ".gz"

        curl_cmd = "curl --output " + gz_fits_file_ref + " \"" + download_url + "\""
        print("curl_cmd =",curl_cmd)

        return_code = os.system(curl_cmd)
        print(f"Command exited with code: {return_code}")

        gunzip_cmd = "gunzip -f " + gz_fits_file_ref
        print("gunzip_cmd =",gunzip_cmd)

        return_code = os.system(gunzip_cmd)
        print(f"Command exited with code: {return_code}")

    except:
        print("*** Error: No FITS filename match found; quitting...")
        exit(64)


    util.compute_image_overlap_area(w_sci,x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,fits_file_ref)


    exit(0)





