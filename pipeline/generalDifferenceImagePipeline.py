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
import math
import numpy as np
from datetime import datetime, timezone, UTC
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
cfg_filename_only = "generalDifferenceImagePipeline.ini"

print("swname =", swname)
print("swvers =", swvers)
print("cfg_filename_only =", cfg_filename_only)


# Compute processing datetime (UT) and processing datetime (Pacific time).

datetime_utc_now = datetime.now(UTC)
proc_utc_datetime = datetime_utc_now.strftime('%Y-%m-%dT%H:%M:%SZ')
datetime_pt_now = datetime_utc_now.replace(tzinfo=timezone.utc).astimezone(tz=to_zone)
proc_pt_datetime_started = datetime_pt_now.strftime('%Y-%m-%dT%H:%M:%S PT')

print("proc_utc_datetime =",proc_utc_datetime)
print("proc_pt_datetime_started =",proc_pt_datetime_started)


# Specify science image.

input_fits_file = 'ADP.2022-07-27T14_56_30.297.fits'
hdr_num = 1                                                    # Second HDU


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

verbose = int(config_input['JOB_PARAMS']['verbose'])
debug = int(config_input['JOB_PARAMS']['debug'])
job_info_s3_bucket_base = config_input['JOB_PARAMS']['job_info_s3_bucket_base']
product_s3_bucket_base = config_input['JOB_PARAMS']['product_s3_bucket_base']
job_config_filename_base = config_input['JOB_PARAMS']['job_config_filename_base']
product_config_filename_base = config_input['JOB_PARAMS']['product_config_filename_base']
refimage_psf_s3_bucket_dir = config_input['JOB_PARAMS']['refimage_psf_s3_bucket_dir']
refimage_psf_filename = config_input['JOB_PARAMS']['refimage_psf_filename']

sca_gain = float(config_input['INSTRUMENT']['sca_gain'])

ppid = int(config_input['SCI_IMAGE']['ppid'])
saturation_level_sciimage = config_input['SCI_IMAGE']['saturation_level']

ppid_refimage = int(config_input['REF_IMAGE']['ppid_refimage'])
max_n_images_to_coadd = int(config_input['REF_IMAGE']['max_n_images_to_coadd'])
naxis1_refimage = int(config_input['REF_IMAGE']['naxis1_refimage'])
naxis2_refimage = int(config_input['REF_IMAGE']['naxis2_refimage'])
cdelt1_refimage = float(config_input['REF_IMAGE']['cdelt1_refimage'])
cdelt2_refimage = float(config_input['REF_IMAGE']['cdelt2_refimage'])
crota2_refimage = float(config_input['REF_IMAGE']['crota2_refimage'])


# Get the awaicgen parameters.  Some of these parameters will be overwritten by this script.
# Do not convert to numerical types, since these will just be passed through (except for those
# overwritten by this script).

awaicgen_dict = {}

awaicgen_dict["awaicgen_input_images_list_file"] = config_input['AWAICGEN']['awaicgen_input_images_list_file']
awaicgen_dict["awaicgen_input_uncert_list_file"] = config_input['AWAICGEN']['awaicgen_input_uncert_list_file']
awaicgen_dict["awaicgen_mosaic_size_x"] = config_input['AWAICGEN']['awaicgen_mosaic_size_x']
awaicgen_dict["awaicgen_mosaic_size_y"] = config_input['AWAICGEN']['awaicgen_mosaic_size_y']
awaicgen_dict["awaicgen_RA_center"] = config_input['AWAICGEN']['awaicgen_RA_center']
awaicgen_dict["awaicgen_Dec_center"] = config_input['AWAICGEN']['awaicgen_Dec_center']
awaicgen_dict["awaicgen_mosaic_rotation"] = config_input['AWAICGEN']['awaicgen_mosaic_rotation']
awaicgen_dict["awaicgen_pixelscale_factor"] = config_input['AWAICGEN']['awaicgen_pixelscale_factor']
awaicgen_dict["awaicgen_pixelscale_absolute"] = config_input['AWAICGEN']['awaicgen_pixelscale_absolute']
awaicgen_dict["awaicgen_mos_cellsize_factor"] = config_input['AWAICGEN']['awaicgen_mos_cellsize_factor']
awaicgen_dict["awaicgen_drizzle_factor"] = config_input['AWAICGEN']['awaicgen_drizzle_factor']
awaicgen_dict["awaicgen_inv_var_weight_flag"] = config_input['AWAICGEN']['awaicgen_inv_var_weight_flag']
awaicgen_dict["awaicgen_pixelflux_scale_flag"] = config_input['AWAICGEN']['awaicgen_pixelflux_scale_flag']
awaicgen_dict["awaicgen_simple_coadd_flag"] = config_input['AWAICGEN']['awaicgen_simple_coadd_flag']
awaicgen_dict["awaicgen_num_threads"] = config_input['AWAICGEN']['awaicgen_num_threads']
awaicgen_dict["awaicgen_unc_sigfigs_retained"] = config_input['AWAICGEN']['awaicgen_unc_sigfigs_retained']
awaicgen_dict["awaicgen_output_mosaic_image_file"] = config_input['AWAICGEN']['awaicgen_output_mosaic_image_file']
awaicgen_dict["awaicgen_output_mosaic_cov_map_file"] = config_input['AWAICGEN']['awaicgen_output_mosaic_cov_map_file']
awaicgen_dict["awaicgen_output_mosaic_uncert_image_file"] = config_input['AWAICGEN']['awaicgen_output_mosaic_uncert_image_file']
awaicgen_dict["awaicgen_debug"] = config_input['AWAICGEN']['awaicgen_debug']
awaicgen_dict["awaicgen_verbose"] = config_input['AWAICGEN']['awaicgen_verbose']


# Update the awaicgen dictionary for quantities that do not vary with sky location.

pixel_scale = math.fabs(cdelt1_refimage)
awaicgen_mosaic_size_x = pixel_scale * float(naxis1_refimage)
awaicgen_mosaic_size_y = pixel_scale * float(naxis2_refimage)

awaicgen_dict["awaicgen_mosaic_size_x"] = str(awaicgen_mosaic_size_x)
awaicgen_dict["awaicgen_mosaic_size_y"] = str(awaicgen_mosaic_size_y)
awaicgen_dict["awaicgen_mosaic_rotation"] = str(crota2_refimage)


# Get the ZOGY parameters.
# Do not convert to numerical types, since these will just be passed through.

zogy_dict = {}

zogy_dict["astrometric_uncert_x"] = config_input['ZOGY']['astrometric_uncert_x']
zogy_dict["astrometric_uncert_y"] = config_input['ZOGY']['astrometric_uncert_y']
zogy_dict["zogy_output_diffimage_file"] = config_input['ZOGY']['zogy_output_diffimage_file']
zogy_dict["post_zogy_keep_diffimg_lower_cov_map_thresh"] = config_input['ZOGY']['post_zogy_keep_diffimg_lower_cov_map_thresh']


# Get the swarp parameters.  Some of these parameters will be overwritten by this script.
# Do not convert to numerical types, since these will just be passed through.

swarp_dict = {}

swarp_dict["swarp_input_image"] = config_input['SWARP']['swarp_input_image']
swarp_dict["swarp_IMAGEOUT_NAME"] = config_input['SWARP']['swarp_IMAGEOUT_NAME']
swarp_dict["swarp_WEIGHTOUT_NAME"] = config_input['SWARP']['swarp_WEIGHTOUT_NAME']
swarp_dict["swarp_HEADER_ONLY"] = config_input['SWARP']['swarp_HEADER_ONLY']
swarp_dict["swarp_HEADER_SUFFIX"] = config_input['SWARP']['swarp_HEADER_SUFFIX']
swarp_dict["swarp_WEIGHT_TYPE"] = config_input['SWARP']['swarp_WEIGHT_TYPE']
swarp_dict["swarp_RESCALE_WEIGHTS"] = config_input['SWARP']['swarp_RESCALE_WEIGHTS']
swarp_dict["swarp_WEIGHT_SUFFIX"] = config_input['SWARP']['swarp_WEIGHT_SUFFIX']
swarp_dict["swarp_WEIGHT_IMAGE"] = config_input['SWARP']['swarp_WEIGHT_IMAGE']
swarp_dict["swarp_WEIGHT_THRESH"] = config_input['SWARP']['swarp_WEIGHT_THRESH']
swarp_dict["swarp_COMBINE"] = config_input['SWARP']['swarp_COMBINE']
swarp_dict["swarp_COMBINE_TYPE"] = config_input['SWARP']['swarp_COMBINE_TYPE']
swarp_dict["swarp_CLIP_AMPFRAC"] = config_input['SWARP']['swarp_CLIP_AMPFRAC']
swarp_dict["swarp_CLIP_SIGMA"] = config_input['SWARP']['swarp_CLIP_SIGMA']
swarp_dict["swarp_CLIP_WRITELOG"] = config_input['SWARP']['swarp_CLIP_WRITELOG']
swarp_dict["swarp_CLIP_LOGNAME"] = config_input['SWARP']['swarp_CLIP_LOGNAME']
swarp_dict["swarp_BLANK_BADPIXELS"] = config_input['SWARP']['swarp_BLANK_BADPIXELS']
swarp_dict["swarp_CELESTIAL_TYPE"] = config_input['SWARP']['swarp_CELESTIAL_TYPE']
swarp_dict["swarp_PROJECTION_TYPE"] = config_input['SWARP']['swarp_PROJECTION_TYPE']
swarp_dict["swarp_PROJECTION_ERR"] = config_input['SWARP']['swarp_PROJECTION_ERR']
swarp_dict["swarp_CENTER_TYPE"] = config_input['SWARP']['swarp_CENTER_TYPE']
swarp_dict["swarp_CENTER"] = config_input['SWARP']['swarp_CENTER']
swarp_dict["swarp_PIXELSCALE_TYPE"] = config_input['SWARP']['swarp_PIXELSCALE_TYPE']
swarp_dict["swarp_PIXEL_SCALE"] = config_input['SWARP']['swarp_PIXEL_SCALE']
swarp_dict["swarp_IMAGE_SIZE"] = config_input['SWARP']['swarp_IMAGE_SIZE']
swarp_dict["swarp_RESAMPLE"] = config_input['SWARP']['swarp_RESAMPLE']
swarp_dict["swarp_RESAMPLE_DIR"] = config_input['SWARP']['swarp_RESAMPLE_DIR']
swarp_dict["swarp_RESAMPLE_SUFFIX"] = config_input['SWARP']['swarp_RESAMPLE_SUFFIX']
swarp_dict["swarp_RESAMPLING_TYPE"] = config_input['SWARP']['swarp_RESAMPLING_TYPE']
swarp_dict["swarp_OVERSAMPLING"] = config_input['SWARP']['swarp_OVERSAMPLING']
swarp_dict["swarp_INTERPOLATE"] = config_input['SWARP']['swarp_INTERPOLATE']
swarp_dict["swarp_FSCALASTRO_TYPE"] = config_input['SWARP']['swarp_FSCALASTRO_TYPE']
swarp_dict["swarp_FSCALE_KEYWORD"] = config_input['SWARP']['swarp_FSCALE_KEYWORD']
swarp_dict["swarp_FSCALE_DEFAULT"] = config_input['SWARP']['swarp_FSCALE_DEFAULT']
swarp_dict["swarp_GAIN_KEYWORD"] = config_input['SWARP']['swarp_GAIN_KEYWORD']
swarp_dict["swarp_GAIN_DEFAULT"] = config_input['SWARP']['swarp_GAIN_DEFAULT']
swarp_dict["swarp_SATLEV_KEYWORD"] = config_input['SWARP']['swarp_SATLEV_KEYWORD']
swarp_dict["swarp_SATLEV_DEFAULT"] = config_input['SWARP']['swarp_SATLEV_DEFAULT']
swarp_dict["swarp_SUBTRACT_BACK"] = config_input['SWARP']['swarp_SUBTRACT_BACK']
swarp_dict["swarp_BACK_TYPE"] = config_input['SWARP']['swarp_BACK_TYPE']
swarp_dict["swarp_BACK_DEFAULT"] = config_input['SWARP']['swarp_BACK_DEFAULT']
swarp_dict["swarp_BACK_SIZE"] = config_input['SWARP']['swarp_BACK_SIZE']
swarp_dict["swarp_BACK_FILTERSIZE"] = config_input['SWARP']['swarp_BACK_FILTERSIZE']
swarp_dict["swarp_BACK_FILTTHRESH"] = config_input['SWARP']['swarp_BACK_FILTTHRESH']
swarp_dict["swarp_VMEM_DIR"] = config_input['SWARP']['swarp_VMEM_DIR']
swarp_dict["swarp_VMEM_MAX"] = config_input['SWARP']['swarp_VMEM_MAX']
swarp_dict["swarp_MEM_MAX"] = config_input['SWARP']['swarp_MEM_MAX']
swarp_dict["swarp_COMBINE_BUFSIZE"] = config_input['SWARP']['swarp_COMBINE_BUFSIZE']
swarp_dict["swarp_DELETE_TMPFILES"] = config_input['SWARP']['swarp_DELETE_TMPFILES']
swarp_dict["swarp_COPY_KEYWORDS"] = config_input['SWARP']['swarp_COPY_KEYWORDS']
swarp_dict["swarp_WRITE_FILEINFO"] = config_input['SWARP']['swarp_WRITE_FILEINFO']
swarp_dict["swarp_WRITE_XML"] = config_input['SWARP']['swarp_WRITE_XML']
swarp_dict["swarp_VERBOSE_TYPE"] = config_input['SWARP']['swarp_VERBOSE_TYPE']
swarp_dict["swarp_NNODES"] = config_input['SWARP']['swarp_NNODES']
swarp_dict["swarp_NODE_INDEX"] = config_input['SWARP']['swarp_NODE_INDEX']
swarp_dict["swarp_NTHREADS"] = config_input['SWARP']['swarp_NTHREADS']
swarp_dict["swarp_NOPENFILES_MAX"] = config_input['SWARP']['swarp_NOPENFILES_MAX']


# Get the sextractor parameters.  Some of these parameters will be overwritten by this script.
# Do not convert to numerical types, since these will just be passed through.

sextractor_diffimage_dict = config_input['SEXTRACTOR_DIFFIMAGE']
sextractor_sciimage_dict = config_input['SEXTRACTOR_SCIIMAGE']

sextractor_refimage_dict = {}
for key in config_input['SEXTRACTOR_REFIMAGE'].keys():
    #print('Input SEXTRACTOR_REFIMAGE: key, value =',key,config_input['SEXTRACTOR_REFIMAGE'][key])
    sextractor_refimage_dict[key] = config_input['SEXTRACTOR_REFIMAGE'][key]

bkgest_dict = config_input['BKGEST']

gainmatch_dict = config_input['GAINMATCH']
psfcat_diffimage_dict = config_input['PSFCAT_DIFFIMAGE']

sextractor_gainmatch_dict = {}
for key in config_input['SEXTRACTOR_GAINMATCH'].keys():
    #print('Input SEXTRACTOR_GAINMATCH: key, value =',key,config_input['SEXTRACTOR_GAINMATCH'][key])
    sextractor_gainmatch_dict[key] = config_input['SEXTRACTOR_GAINMATCH'][key]

sfft_dict = config_input['SFFT']

naive_diffimage_dict = config_input['NAIVE_DIFFIMAGE']


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


    util.compute_image_overlap_area(w_sci,naxis1,naxis2,x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,fits_file_ref)


    exit(0)





