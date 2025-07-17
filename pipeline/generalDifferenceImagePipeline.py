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


    # Query for 2MASS images that overlap sky position within 1.0 arcseconds.

    twomass_service = vo.dal.SIAService("https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?type=at&ds=asky&")
    im_table = twomass_service.search(pos=pos, size=1.0*u.arcsec)
    print(im_table.to_table())


    # Columns in table.
    # name,download,center_ra,center_dec,naxes,naxis,scale,format,crpix,crval,crota2,band,bref,bhi,blo,pers_art,glint_art,type,dataset,pixflags,id,scntr,date,hem,scan,image,ut_date,coadd_key,seesh,magzp,msnr10,bin


    # Get the returned 2MASS image(s) for the filter of interest.
    # Pick the one that overlaps the science image the most.

    max_percent_overlap_area = 0.0
    intersecting_polygon_file = "intersecting_polygon.txt"

    for i in range(len(im_table)):

        if "fits" not in im_table[i]['download']:
            continue

        if im_table[i]['band'] != filter_science:
            continue

        print(im_table[i].getdataurl())


        # Download the 2MASS image and decompress it.
        #
        # curl --output hi0600256.fits.gz "https://irsa.ipac.caltech.edu:443/cgi-bin/2MASS/IM/nph-im?ds=asky&atdir=/ti08&dh=000616n&scan=060&name=hi0600256.fits"
        # gunzip hi0600256.fits.gz

        download_url = im_table[i].getdataurl()

        filename_match = re.match(r".+?name\=(.+?.fits)", download_url)

        try:
            provisional_fits_file_ref = filename_match.group(1)
            print("provisional_fits_file_ref =",provisional_fits_file_ref)

            gz_provisional_fits_file_ref = provisional_fits_file_ref + ".gz"

            curl_cmd = "curl --output " + gz_provisional_fits_file_ref + " \"" + download_url + "\""
            print("curl_cmd =",curl_cmd)

            return_code = os.system(curl_cmd)
            print(f"Command exited with code: {return_code}")

            gunzip_cmd = "gunzip -f " + gz_provisional_fits_file_ref
            print("gunzip_cmd =",gunzip_cmd)

            return_code = os.system(gunzip_cmd)
            print(f"Command exited with code: {return_code}")

        except:
            print("*** Error: No FITS filename match found; quitting...")
            exit(64)


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


    with fits.open(fits_file_ref) as hdul:

        nframes_refimage = hdul[0].header["NUMFRMS"]



    # Get job configuration parameters.

    verbose = int(config_input['JOB_PARAMS']['verbose'])
    debug = int(config_input['JOB_PARAMS']['debug'])
    refimage_psf_filename = config_input['JOB_PARAMS']['refimage_psf_filename']

    product_config_filename_base = config_input['JOB_PARAMS']['product_config_filename_base']

    sca_gain = float(config_input['INSTRUMENT']['sca_gain'])

    ppid_sciimage = int(config_input['SCI_IMAGE']['ppid'])
    saturation_level_sciimage = float(config_input['SCI_IMAGE']['saturation_level'])


    gain_refimage = float(config_input['REF_IMAGE']['gain_refimage'])
    ppid_refimage = int(config_input['REF_IMAGE']['ppid_refimage'])
    max_n_images_to_coadd = int(config_input['REF_IMAGE']['max_n_images_to_coadd'])
    naxis1_refimage = int(config_input['REF_IMAGE']['naxis1_refimage'])
    naxis2_refimage = int(config_input['REF_IMAGE']['naxis2_refimage'])
    cdelt1_refimage = float(config_input['REF_IMAGE']['cdelt1_refimage'])
    cdelt2_refimage = float(config_input['REF_IMAGE']['cdelt2_refimage'])
    crota2_refimage = float(config_input['REF_IMAGE']['crota2_refimage'])

    astrometric_uncert_x = float(config_input['ZOGY']['astrometric_uncert_x'])
    astrometric_uncert_y = float(config_input['ZOGY']['astrometric_uncert_y'])
    zogy_output_diffimage_file = config_input['ZOGY']['zogy_output_diffimage_file']
    post_zogy_keep_diffimg_lower_cov_map_thresh = float(config_input['ZOGY']['post_zogy_keep_diffimg_lower_cov_map_thresh'])

    awaicgen_dict = config_input['AWAICGEN']

    swarp_dict = config_input['SWARP']

    sextractor_diffimage_dict = config_input['SEXTRACTOR_DIFFIMAGE']
    sextractor_sciimage_dict = config_input['SEXTRACTOR_SCIIMAGE']
    sextractor_refimage_dict = config_input['SEXTRACTOR_REFIMAGE']
    bkgest_dict = config_input['BKGEST']
    gainmatch_dict = config_input['GAINMATCH']
    psfcat_diffimage_dict = config_input['PSFCAT_DIFFIMAGE']
    sextractor_gainmatch_dict = config_input['SEXTRACTOR_GAINMATCH']
    sfft_dict = config_input['SFFT']
    naive_diffimage_dict = config_input['NAIVE_DIFFIMAGE']

    print("max_n_images_to_coadd =", max_n_images_to_coadd)

    saturation_level_refimage = float(sextractor_refimage_dict["sextractor_SATUR_LEVEL".lower()])



    # Code-timing benchmark.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds after downloading reference image =",end_time_benchmark - start_time_benchmark)
    start_time_benchmark = end_time_benchmark



    # Compute reference-image statistics.

    n_sigma = 3.0
    hdu_index = 0

    stats_refimage = util.fits_data_statistics_with_clipping(fits_file_ref,\
                                                             n_sigma,\
                                                             hdu_index,\
                                                             saturation_level_refimage)

    avg_refimage = stats_refimage["clippedavg"]
    std_refimage = stats_refimage["clippedstd"]
    cnt_refimage = stats_refimage["nkept"]
    noutliers_refimage = stats_refimage["noutliers"]
    gmed_refimage = stats_refimage["gmed"]
    datascale_refimage = stats_refimage["gsigma"]
    gmin_refimage = stats_refimage["gdatamin"]
    gmax_refimage = stats_refimage["gdatamax"]
    npixsat_refimage = stats_refimage["satcount"]
    npixnan_refimage = stats_refimage["nancount"]


    # Compute reference-image uncertainty image via simple model (photon noise only).

    fits_file_ref_unc = fits_file_ref.replace(".fits","_unc.fits")

    util.compute_uncertainty_image_via_simple_model(fits_file_ref,
                                                    hdu_index,
                                                    fits_file_ref_unc,
                                                    gain_refimage,
                                                    nframes_refimage)


    # Compute reference-image coverage map.

    fits_file_ref_cov = fits_file_ref.replace(".fits","_cov.fits")

    rfis.compute_coverage_map(fits_file_ref,
                              hdu_index,
                              fits_file_ref_cov,
                              nframes_refimage)


    # Generate reference-image catalog.

    filename_refimage_catalog = fits_file_ref.replace(".fits","_secat.txt")

    rfis.generateReferenceImageCatalog(fits_file_ref,
                                       fits_file_ref_unc,
                                       cfg_path,
                                       sextractor_refimage_dict,
                                       filename_refimage_catalog)


    # Compute additional quantities needed for the reference-image PSF.

    sextractor_refimage_paramsfile = cfg_path + "/srcExtractParamsRefImage.inp";
    params_to_get_refimage = ["FWHM_IMAGE"]

    vals_refimage = util.parse_ascii_text_sextractor_catalog(filename_refimage_catalog,
                                                             sextractor_refimage_paramsfile,
                                                             params_to_get_refimage)

    nsexcatsources_refimage = len(vals_refimage)

    fwhm_ref_vals = []
    for val in vals_refimage:
        fwhm_ref_vals.append(float(val[0]))

    np_fwhm_ref_vals = np.array(fwhm_ref_vals)

    fwhm_ref_minpix = np.nanmin(np_fwhm_ref_vals)
    fwhm_ref_maxpix = np.nanmax(np_fwhm_ref_vals)
    fwhm_ref_medpix = np.nanmedian(np_fwhm_ref_vals)

    print("fwhm_ref_medpix,fwhm_ref_minpix,fwhm_ref_maxpix =",fwhm_ref_medpix,fwhm_ref_minpix,fwhm_ref_maxpix)

    fwhm_ref = fwhm_ref_medpix
    if fwhm_ref < 0.0:
        fwhm_ref = 2.0

    print("fwhm_ref =",fwhm_ref)


    # Code-timing benchmark.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds after downloading or generating reference image =",
        end_time_benchmark - start_time_benchmark)
    start_time_benchmark = end_time_benchmark


    # Populate config-file dictionary for products.

    product_config_ini_filename = product_config_filename_base + str(jid) + ".ini"

    product_config = configparser.ConfigParser()

    product_config['JOB_PARAMS'] = {'debug': str(debug),
                                 'swname': swname,
                                 'swvers': swvers}

    product_config['JOB_PARAMS']['jid'] = str(jid)
    product_config['JOB_PARAMS']['verbose'] = str(verbose)
    product_config['JOB_PARAMS']['job_started'] = str(proc_pt_datetime_started)


    product_config['REF_IMAGE'] = {}

    product_config['REF_IMAGE']['ppid'] = str(ppid_sciimage)


    product_config['REF_IMAGE']['nframes'] = str(nframes_refimage)
    product_config['REF_IMAGE']['npixsat'] = str(npixsat_refimage)
    product_config['REF_IMAGE']['npixnan'] = str(npixnan_refimage)
    product_config['REF_IMAGE']['clmean'] = str(avg_refimage)
    product_config['REF_IMAGE']['clstddev'] = str(std_refimage)
    product_config['REF_IMAGE']['clnoutliers'] = str(noutliers_refimage)
    product_config['REF_IMAGE']['gmedian'] = str(gmed_refimage)
    product_config['REF_IMAGE']['datascale'] = str(datascale_refimage)
    product_config['REF_IMAGE']['gmin'] = str(gmin_refimage)
    product_config['REF_IMAGE']['gmax'] = str(gmax_refimage)
    product_config['REF_IMAGE']['fwhm_ref_medpix'] = str(fwhm_ref_medpix)
    product_config['REF_IMAGE']['fwhm_ref_minpix'] = str(fwhm_ref_minpix)
    product_config['REF_IMAGE']['fwhm_ref_maxpix'] = str(fwhm_ref_maxpix)
    product_config['REF_IMAGE']['nsexcatsources'] = str(nsexcatsources_refimage)


    # Compute image statistics for image resizing.

    n_sigma = 3.0

    stats_sci_img = util.fits_data_statistics_with_clipping(filename_science_image,
                                                            n_sigma,
                                                            hdu_index_science,
                                                            saturation_level_sciimage)

    avg_sci_img = stats_sci_img["clippedavg"]


    # Reformat the science image FITS file
    # so that the image data are contained in the PRIMARY header.
    # Compute uncertainty image via simple model (photon noise only).
    # Resize images
    # from
    # NAXIS1  =                 2099 / length of data axis 1
    # NAXIS2  =                 2100 / length of data axis 2
    # to
    # NAXIS1  =                 2099 / length of data axis 1
    # NAXIS2  =                 2101 / length of data axis 2
    # (odd number of pixels on each side).

    reformatted_science_image_filename = filename_science_image.replace(".fits","_reformatted.fits")
    reformatted_science_uncert_image_filename = filename_science_image.replace(".fits","_reformatted_unc.fits")

    append_extra_col = True
    num_extra_cols = 2
    append_extra_row = True
    num_extra_rows = 1

    dfis.reformat_science_fits_file_and_compute_uncertainty_image_via_simple_model(filename_science_image,
                                                                                   hdu_index_science,
                                                                                   reformatted_science_image_filename,
                                                                                   reformatted_science_uncert_image_filename,
                                                                                   append_extra_col,
                                                                                   num_extra_cols,
                                                                                   append_extra_row,
                                                                                   num_extra_rows,
                                                                                   sca_gain,
                                                                                   avg_sci_img)








    # Generate science-image catalog.

    filename_sciimage_catalog = reformatted_science_image_filename.replace(".fits","_secat.txt")

    util.generateScienceImageCatalog(reformatted_science_image_filename,
                                                                                   reformatted_science_uncert_image_filename,
                                                                                   cfg_path,
                                                                                   sextractor_sciimage_dict,
                                                                                   filename_sciimage_catalog)


    # Compute additional quantities needed for the science-image PSF.

    sextractor_sciimage_paramsfile = cfg_path + "/srcExtractParamsSciImage.inp";
    params_to_get_sciimage = ["FWHM_IMAGE"]

    vals_sciimage = util.parse_ascii_text_sextractor_catalog(filename_sciimage_catalog,
                                                             sextractor_sciimage_paramsfile,
                                                             params_to_get_sciimage)

    nsexcatsources_sciimage = len(vals_sciimage)

    fwhm_sci_vals = []
    for val in vals_sciimage:
        fwhm_sci_vals.append(float(val[0]))

    np_fwhm_sci_vals = np.array(fwhm_sci_vals)

    fwhm_sci_minpix = np.nanmin(np_fwhm_sci_vals)
    fwhm_sci_maxpix = np.nanmax(np_fwhm_sci_vals)
    fwhm_sci_medpix = np.nanmedian(np_fwhm_sci_vals)

    print("fwhm_sci_medpix,fwhm_sci_minpix,fwhm_sci_maxpix =",fwhm_sci_medpix,fwhm_sci_minpix,fwhm_sci_maxpix)

    fwhm_sci = fwhm_sci_medpix
    if fwhm_sci < 0.0:
        fwhm_sci = 2.0

    print("fwhm_sci =",fwhm_sci)









    # Swarp the reference image and associated uncertainty image into the distortion frame of the science image.
    # Since the reference image was made by awaicgen, there is no geometric image distortion,
    # and, hence, no need to convert from sip to pv distortion, so the following flag is set to False.
    # Set the following flag to True only for the case where the reference image is a single Roman SCA image.

    hdu_index_for_science_image_data = 0
    hdu_index_for_reference_image_data = 0
    pv_convert_flag_for_science_image_data = False                   # TODO
    pv_convert_flag_for_reference_image_data = False                   # TODO

    sci_fits_file_with_pv,\
        ref_fits_file_with_pv,\
        ref_cov_fits_file_with_pv,\
        ref_uncert_fits_file_with_pv,\
        output_resampled_reference_image,\
        output_resampled_reference_cov_map,\
        output_resampled_reference_uncert_image =\
        util.resample_reference_image_to_science_image_with_pv_distortion(reformatted_science_image_filename,\
                                                                          hdu_index_for_science_image_data,\
                                                                          fits_file_ref,\
                                                                          fits_file_ref_cov,\
                                                                          fits_file_ref_unc,\
                                                                          hdu_index_for_reference_image_data,\
                                                                          pv_convert_flag_for_science_image_data,\
                                                                          pv_convert_flag_for_reference_image_data,\
                                                                          swarp_dict)


    # Code-timing benchmark.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds after swarping images =",end_time_benchmark - start_time_benchmark)
    start_time_benchmark = end_time_benchmark


    # Compute image statistics for ZOGY.

    n_sigma = 3.0
    hdu_index = 0

    stats_sci_img = util.fits_data_statistics_with_clipping(reformatted_science_image_filename,\
                                                            n_sigma,\
                                                            hdu_index,\
                                                            saturation_level_sciimage)

    avg_sci_img = stats_sci_img["clippedavg"]
    std_sci_img = stats_sci_img["clippedstd"]
    cnt_sci_img = stats_sci_img["nkept"]

    stats_ref_img = util.fits_data_statistics_with_clipping(output_resampled_reference_image,\
                                                            n_sigma,\
                                                            hdu_index,\
                                                            saturation_level_refimage)

    avg_ref_img = stats_ref_img["clippedavg"]
    std_ref_img = stats_ref_img["clippedstd"]
    cnt_ref_img = stats_ref_img["nkept"]

    print("avg_sci_img,std_sci_img,cnt_sci_img =",avg_sci_img,std_sci_img,cnt_sci_img)
    print("avg_ref_img,std_ref_img,cnt_ref_img =",avg_ref_img,std_ref_img,cnt_ref_img)




    # Subtract background from science image.  Since the reference image has been swarped,
    # it already has the background subtracted.

    bkgest_code = rapid_sw + '/c/bin/bkgest'
    bkgest_include_dir = rapid_sw + '/c/include'
    filename_bkg_subbed_science_image = 'bkg_subbed_science_image.fits'
    filename_global_clippedmean_sciimage_tbl = 'global_clippedmean_science_image.tbl'

    bkgest_cmd = [bkgest_code,
                  '-i',
                  reformatted_science_image_filename,
                  '-f',
                  bkgest_dict["output_image_type"],
                  '-c',
                  bkgest_dict["clippedmean_calc_type"],
                  '-g',
                  bkgest_dict["local_clippedmean_grid_spacing"],
                  '-w',
                  bkgest_dict["local_clippedmean_input_window"],
                  '-a',
                  bkgest_include_dir,
                  '-ot',
                  filename_global_clippedmean_sciimage_tbl,
                  '-o2',
                  filename_bkg_subbed_science_image]

    exitcode_from_bkgest = util.execute_command(bkgest_cmd)


    # Code-timing benchmark.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds after running bkgest on science image =",
        end_time_benchmark - start_time_benchmark)
    start_time_benchmark = end_time_benchmark


    ######################################################################################
    # Gain-match science and reference images by generating SExtractor catalogs for each
    # and then computing scale factor.  To apply, multiply reference image by scalefacref.
    ######################################################################################

    scalefac,dxrmsfin,dyrmsfin,dxmedianfin,dymedianfin = dfis.gainMatchScienceAndReferenceImages(filename_bkg_subbed_science_image,
                                                                                                 reformatted_science_uncert_image_filename,
                                                                                                 output_resampled_reference_image,
                                                                                                 output_resampled_reference_uncert_image,
                                                                                                 cfg_path,
                                                                                                 gainmatch_dict,
                                                                                                 sextractor_gainmatch_dict,
                                                                                                 fwhm_sci,
                                                                                                 fwhm_ref,
                                                                                                 astrometric_uncert_x,
                                                                                                 astrometric_uncert_y)

    print("scalefac,dxrmsfin,dyrmsfin,dxmedianfin,dymedianfin =",scalefac,dxrmsfin,dyrmsfin,dxmedianfin,dymedianfin)

    scalefacref = 1. / scalefac


    # Compute resampled gain-matched reference image.

    output_resampled_gainmatched_reference_image = output_resampled_reference_image.replace(".fits","_gainmatched.fits")
    util.scale_image_data(output_resampled_reference_image,scalefacref,output_resampled_gainmatched_reference_image)


    # Code-timing benchmark.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds after running gainMatchScienceAndReferenceImages =",
        end_time_benchmark - start_time_benchmark)
    start_time_benchmark = end_time_benchmark


    # Replace NaNs, if any, in ZOGY input images.  Use the same saturation level rate since they are gain-matched.

    print("saturation_level_sciimage =",saturation_level_sciimage)
    nan_indices_sciimage = util.replace_nans_with_sat_val_rate(filename_bkg_subbed_science_image,saturation_level_sciimage)
    nan_indices_refimage = util.replace_nans_with_sat_val_rate(output_resampled_gainmatched_reference_image,saturation_level_sciimage)


    # Apply subpixel orthogonal offsets to ZOGY input reference image.

    util.apply_subpixel_orthogonal_offsets(output_resampled_gainmatched_reference_image,dxmedianfin,dymedianfin)


    # Code-timing benchmark.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds after replacing NaNs, applying image offsets, etc. =",
        end_time_benchmark - start_time_benchmark)
    start_time_benchmark = end_time_benchmark


    # Generate PSFs for science and reference images.

    nside = 61

    filename_sci_psf = "science_image_psf.fits"
    util.generate_2d_gaussian_psf(fwhm_sci,nside,filename_sci_psf)

    filename_ref_psf = "reference_image_psf.fits"
    util.generate_2d_gaussian_psf(fwhm_ref,nside,filename_ref_psf)


    #################################################################################################################
    # The image data in filename_science_image and sci_fits_file_with_pv FITS files are the same, only the
    # representation of geometric distortion in the FITS headers are different (sip versus pv).
    #
    # ZOGY only cares about the image data, not what is in the FITS headers.
    # Usage: python py_zogy.py <NewImage> <RefImage> <NewPSF> <RefPSF> <NewSigmaImage> <RefSigmaImage>
    #                    <NewSigmaMode> <RefSigmaMode> <AstUncertX> <AstUncertY> <DiffImage> <DiffPSF> <ScorrImage>
    #
    # Assume top-level directory of RAPID-pipeline git repo is mapped to rapid_sw inside Docker container.
    #################################################################################################################


    python_cmd = '/usr/bin/python3.12'
    zogy_code = rapid_sw + '/modules/zogy/v21Aug2018/py_zogy.py'
    filename_diffimage = 'diffimage.fits'
    filename_diffpsf = 'diffpsf.fits'
    filename_scorrimage = 'scorrimage.fits'

    zogy_cmd = [python_cmd,
                zogy_code,
                filename_bkg_subbed_science_image,
                output_resampled_gainmatched_reference_image,
                filename_sci_psf,
                filename_ref_psf,
                reformatted_science_uncert_image_filename,
                output_resampled_reference_uncert_image,
                str(std_sci_img),
                str(std_ref_img),
                str(dxrmsfin),
                str(dyrmsfin),
                filename_diffimage,
                filename_diffpsf,
                filename_scorrimage]

    exitcode_from_zogy = util.execute_command(zogy_cmd)


    # Code-timing benchmark.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds after running ZOGY =",
        end_time_benchmark - start_time_benchmark)
    start_time_benchmark = end_time_benchmark


    # Mask difference image with output_resampled_reference_cov_map.

    filename_diffimage_masked = zogy_output_diffimage_file                     # Nominally diffimage_masked.fits
    filename_scorrimage_masked = 'scorrimage_masked.fits'

    dfis.mask_difference_image_with_resampled_reference_cov_map(filename_diffimage,
                                                                output_resampled_reference_cov_map,
                                                                filename_diffimage_masked,
                                                                post_zogy_keep_diffimg_lower_cov_map_thresh)

    dfis.mask_difference_image_with_resampled_reference_cov_map(filename_scorrimage,
                                                                output_resampled_reference_cov_map,
                                                                filename_scorrimage_masked,
                                                                post_zogy_keep_diffimg_lower_cov_map_thresh)


    # Restore NaNs that were masked prior to executing ZOGY.

    if nan_indices_sciimage:
        util.restore_nans(filename_diffimage_masked,nan_indices_sciimage)

    if nan_indices_refimage:
        util.restore_nans(filename_diffimage_masked,nan_indices_refimage)

    if nan_indices_sciimage:
        util.restore_nans(filename_scorrimage_masked,nan_indices_sciimage)

    if nan_indices_refimage:
        util.restore_nans(filename_scorrimage_masked,nan_indices_refimage)


    # Code-timing benchmark.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds after masking ZOGY difference image =",
        end_time_benchmark - start_time_benchmark)
    start_time_benchmark = end_time_benchmark


    # Generate diffimage uncertainty image, which will be the weight image for sextractor_WEIGHT_IMAGE.

    filename_diffimage_unc_masked = 'diffimage_uncert_masked.fits'
    dfis.compute_diffimage_uncertainty(sca_gain,
                                       reformatted_science_image_filename,
                                       output_resampled_gainmatched_reference_image,
                                       output_resampled_reference_cov_map,
                                       filename_diffimage_masked,
                                       filename_diffimage_unc_masked)
    filename_weight_image = filename_diffimage_unc_masked
    filename_diffimage_sextractor_catalog = filename_diffimage_masked.replace(".fits",".txt")


    # Compute SExtractor catalog for masked difference image.
    # Execute SExtractor to first detect candidates on Scorr (S/N) match-filter
    # image, then use to perform aperture phot on difference image to generate
    # raw ascii catalog file.

    sextractor_diffimage_paramsfile = cfg_path + "/srcExtractParamsDiffImage.inp";

    sextractor_diffimage_dict["sextractor_detection_image".lower()] = filename_scorrimage_masked
    sextractor_diffimage_dict["sextractor_input_image".lower()] = filename_diffimage_masked
    sextractor_diffimage_dict["sextractor_WEIGHT_IMAGE".lower()] = filename_weight_image
    sextractor_diffimage_dict["sextractor_PARAMETERS_NAME".lower()] = sextractor_diffimage_paramsfile
    sextractor_diffimage_dict["sextractor_FILTER_NAME".lower()] = cfg_path + "/srcExtractDiffImageFilter.conv"
    sextractor_diffimage_dict["sextractor_STARNNW_NAME".lower()] = cfg_path + "/srcExtractDiffImageStarGalaxyClassifier.nnw"
    sextractor_diffimage_dict["sextractor_CATALOG_NAME".lower()] = filename_diffimage_sextractor_catalog
    sextractor_cmd = util.build_sextractor_command_line_args(sextractor_diffimage_dict)
    exitcode_from_sextractor = util.execute_command(sextractor_cmd)

    params_to_get_diffimage = ["XWIN_IMAGE","YWIN_IMAGE","FLUX_APER_6"]

    vals_diffimage = util.parse_ascii_text_sextractor_catalog(filename_diffimage_sextractor_catalog,
                                                              sextractor_diffimage_paramsfile,
                                                              params_to_get_diffimage)

    nsexcatsources_diffimage = len(vals_diffimage)

    print("nsexcatsources_diffimage =",nsexcatsources_diffimage)


    # Code-timing benchmark.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds after running SExtractor on ZOGY difference image =",
        end_time_benchmark - start_time_benchmark)
    start_time_benchmark = end_time_benchmark







    # Define ZOGY dictionary in config-file dictionary for products.

    product_config['ZOGY'] = {}

    product_config['ZOGY']['ppid'] = str(ppid_sciimage)



    # The following sky positions are correct for the difference image
    # only because the current code reprojects the reference image
    # into the distorted grid of the science image.

    product_config['ZOGY']['ra0'] = str(ra0)
    product_config['ZOGY']['dec0'] = str(dec0)
    product_config['ZOGY']['ra1'] = str(ra1)
    product_config['ZOGY']['dec1'] = str(dec1)
    product_config['ZOGY']['ra2'] = str(ra2)
    product_config['ZOGY']['dec2'] = str(dec2)
    product_config['ZOGY']['ra3'] = str(ra3)
    product_config['ZOGY']['dec3'] = str(dec3)
    product_config['ZOGY']['ra4'] = str(ra4)
    product_config['ZOGY']['dec4'] = str(dec4)

    product_config['ZOGY']['nsexcatsources'] = str(nsexcatsources_diffimage)
    product_config['ZOGY']['scalefacref'] = str(scalefacref)
    product_config['ZOGY']['dxrmsfin'] = str(dxrmsfin)
    product_config['ZOGY']['dyrmsfin'] = str(dyrmsfin)
    product_config['ZOGY']['dxmedianfin'] = str(dxmedianfin)
    product_config['ZOGY']['dymedianfin'] = str(dymedianfin)



    #################################################################################################################
    # Optionally run SFFT to generate an alternate difference image and catalog.
    # Filenames for the SExtractor segmented maps are provided, and, if they do not exist, they will be generated.
    # Output files (constructed by the script, but not provided as input:
    #    sfftdiffimage_masked.fits
    #    sfftsoln.fits
    #################################################################################################################

    run_sfft = eval(sfft_dict['run_sfft'])

    # Always leave as True, and can only be reset to False if and only if SFFT runs and fails.
    run_sfft_was_successful = True

    if run_sfft:

        # Cannot run under python3.11 because scikit-learn fails to install.
        python_cmd = '/usr/bin/python3'
        sfft_code = rapid_sw + '/modules/sfft/sfft_rapid.py'
        filename_scifile = filename_bkg_subbed_science_image
        filename_reffile = output_resampled_gainmatched_reference_image
        filename_scisegm = 'sfftscisegm.fits'
        filename_refsegm = 'sfftrefsegm.fits'

        crossconv_flag = eval(sfft_dict['crossconv_flag'])

        if crossconv_flag:
            filename_sfftdiffimage = 'sfftdiffimage_cconv_masked.fits'
            filename_sfftsoln = 'sfftsoln_cconv.fits'
        else:
            filename_sfftdiffimage = 'sfftdiffimage_masked.fits'
            filename_sfftsoln = 'sfftsoln.fits'

        filename_dcdiff = 'sfftdiffimage_dconv_masked.fits'

        # A quirk in the software requires prepended "./" to input filenames.
        sfft_cmd = [python_cmd,
                    sfft_code,
                    "./" + filename_scifile,
                    "./" + filename_reffile,
                    filename_scisegm,
                    filename_refsegm]

        if crossconv_flag:
            sfft_cmd.append("--crossconv")
            sfft_cmd.append("--scipsf")
            sfft_cmd.append(filename_sci_psf)
            sfft_cmd.append("--refpsf")
            sfft_cmd.append(filename_ref_psf)

        exitcode_from_sfft = util.execute_command(sfft_cmd)

        if int(exitcode_from_sfft) != 0:
            run_sfft_was_successful = False


        # Code-timing benchmark.

        end_time_benchmark = time.time()
        print("Elapsed time in seconds after running SFFT =",
            end_time_benchmark - start_time_benchmark)
        start_time_benchmark = end_time_benchmark

        if run_sfft_was_successful:

            # Generate SFFT diffimage uncertainty image, which will be the weight image for sextractor_WEIGHT_IMAGE.

            filename_sfftdiffimage_unc = 'sfftdiffimage_uncert_masked.fits'
            dfis.compute_diffimage_uncertainty(sca_gain,
                                               reformatted_science_image_filename,
                                               output_resampled_gainmatched_reference_image,
                                               output_resampled_reference_cov_map,
                                               filename_sfftdiffimage,
                                               filename_sfftdiffimage_unc)
            filename_weight_image = filename_sfftdiffimage_unc
            filename_sfftdiffimage_sextractor_catalog = filename_sfftdiffimage.replace(".fits",".txt")


            # Compute SExtractor catalog for masked difference image.
            # Execute SExtractor to first detect candidates on Scorr (S/N) match-filter
            # image, then use to perform aperture phot on difference image to generate
            # raw ascii catalog file.

            sextractor_diffimage_paramsfile = cfg_path + "/srcExtractParamsDiffImage.inp";

            sextractor_diffimage_dict["sextractor_detection_image".lower()] = filename_sfftdiffimage
            sextractor_diffimage_dict["sextractor_input_image".lower()] = filename_sfftdiffimage
            sextractor_diffimage_dict["sextractor_WEIGHT_IMAGE".lower()] = filename_weight_image
            sextractor_diffimage_dict["sextractor_PARAMETERS_NAME".lower()] = sextractor_diffimage_paramsfile
            sextractor_diffimage_dict["sextractor_FILTER_NAME".lower()] = cfg_path + "/srcExtractDiffImageFilter.conv"
            sextractor_diffimage_dict["sextractor_STARNNW_NAME".lower()] = cfg_path + "/srcExtractDiffImageStarGalaxyClassifier.nnw"
            sextractor_diffimage_dict["sextractor_CATALOG_NAME".lower()] = filename_sfftdiffimage_sextractor_catalog
            sextractor_cmd = util.build_sextractor_command_line_args(sextractor_diffimage_dict)
            exitcode_from_sextractor = util.execute_command(sextractor_cmd)

            params_to_get_diffimage = ["XWIN_IMAGE","YWIN_IMAGE","FLUX_APER_6"]

            vals_sfftdiffimage = util.parse_ascii_text_sextractor_catalog(filename_sfftdiffimage_sextractor_catalog,
                                                                          sextractor_diffimage_paramsfile,
                                                                          params_to_get_diffimage)

            nsexcatsources_sfftdiffimage = len(vals_sfftdiffimage)

            print("nsexcatsources_sfftdiffimage =",nsexcatsources_sfftdiffimage)


            # Code-timing benchmark.

            end_time_benchmark = time.time()
            print("Elapsed time in seconds after running SExtractor on SFFT difference image =",
                end_time_benchmark - start_time_benchmark)
            start_time_benchmark = end_time_benchmark




    #################################################################################################################
    # Compute naive image difference.
    #################################################################################################################

    naive_diffimage_flag = eval(naive_diffimage_dict['naive_diffimage_flag'])

    if naive_diffimage_flag:

        filename_naive_diffimage = "naive_diffimage.fits"

        util.compute_naive_difference_image(filename_bkg_subbed_science_image,
                                            output_resampled_gainmatched_reference_image,
                                            filename_naive_diffimage)


        # Mask naive difference image with output_resampled_reference_cov_map.

        filename_naive_diffimage_masked = naive_diffimage_dict['naive_output_diffimage_file']

        dfis.mask_difference_image_with_resampled_reference_cov_map(filename_naive_diffimage,
                                                                    output_resampled_reference_cov_map,
                                                                    filename_naive_diffimage_masked,
                                                                    post_zogy_keep_diffimg_lower_cov_map_thresh)


        # Restore NaNs that were masked prior to executing ZOGY.

        if nan_indices_sciimage:
            util.restore_nans(filename_naive_diffimage_masked,nan_indices_sciimage)

        if nan_indices_refimage:
            util.restore_nans(filename_naive_diffimage_masked,nan_indices_refimage)


        # Code-timing benchmark.

        end_time_benchmark = time.time()
        print("Elapsed time in seconds after computing naive image difference =",
            end_time_benchmark - start_time_benchmark)
        start_time_benchmark = end_time_benchmark


        # Generate naive diffimage uncertainty image, which will be the weight image for sextractor_WEIGHT_IMAGE.

        filename_naivediffimage_unc = filename_naive_diffimage_masked.replace("masked.fits","uncert_masked.fits")

        dfis.compute_diffimage_uncertainty(sca_gain,
                                           reformatted_science_image_filename,
                                           output_resampled_gainmatched_reference_image,
                                           output_resampled_reference_cov_map,
                                           filename_naive_diffimage_masked,
                                           filename_naivediffimage_unc)
        filename_weight_image = filename_naivediffimage_unc
        filename_naivediffimage_sextractor_catalog = filename_naive_diffimage_masked.replace(".fits",".txt")


        # Compute SExtractor catalog for masked difference image.
        # Execute SExtractor to first detect candidates on Scorr (S/N) match-filter
        # image, then use to perform aperture phot on difference image to generate
        # raw ascii catalog file.

        sextractor_diffimage_paramsfile = cfg_path + "/srcExtractParamsDiffImage.inp";

        sextractor_diffimage_dict["sextractor_detection_image".lower()] = filename_naive_diffimage_masked
        sextractor_diffimage_dict["sextractor_input_image".lower()] = filename_naive_diffimage_masked
        sextractor_diffimage_dict["sextractor_WEIGHT_IMAGE".lower()] = filename_weight_image
        sextractor_diffimage_dict["sextractor_PARAMETERS_NAME".lower()] = sextractor_diffimage_paramsfile
        sextractor_diffimage_dict["sextractor_FILTER_NAME".lower()] = cfg_path + "/srcExtractDiffImageFilter.conv"
        sextractor_diffimage_dict["sextractor_STARNNW_NAME".lower()] = cfg_path + "/srcExtractDiffImageStarGalaxyClassifier.nnw"
        sextractor_diffimage_dict["sextractor_CATALOG_NAME".lower()] = filename_naivediffimage_sextractor_catalog
        sextractor_cmd = util.build_sextractor_command_line_args(sextractor_diffimage_dict)
        exitcode_from_sextractor = util.execute_command(sextractor_cmd)

        params_to_get_diffimage = ["XWIN_IMAGE","YWIN_IMAGE","FLUX_APER_6"]

        vals_naivediffimage = util.parse_ascii_text_sextractor_catalog(filename_naivediffimage_sextractor_catalog,
                                                                      sextractor_diffimage_paramsfile,
                                                                      params_to_get_diffimage)

        nsexcatsources_naivediffimage = len(vals_naivediffimage)

        print("nsexcatsources_naivediffimage =",nsexcatsources_naivediffimage)


        # Code-timing benchmark.

        end_time_benchmark = time.time()
        print("Elapsed time in seconds after running SExtractor on SFFT difference image =",
            end_time_benchmark - start_time_benchmark)
        start_time_benchmark = end_time_benchmark


    # Get listing of working directory as a diagnostic.

    ls_cmd = ['ls','-ltr']
    exitcode_from_ls = util.execute_command(ls_cmd)


    # Get timestamp job ended in Pacific Time for Jobs database record later.

    datetime_utc_now = datetime.now(UTC)
    proc_utc_datetime = datetime_utc_now.strftime('%Y-%m-%dT%H:%M:%SZ')
    datetime_pt_now = datetime_utc_now.replace(tzinfo=timezone.utc).astimezone(tz=to_zone)
    proc_pt_datetime_ended = datetime_pt_now.strftime('%Y-%m-%dT%H:%M:%S PT')

    print("proc_pt_datetime_ended =",proc_pt_datetime_ended)

    product_config['JOB_PARAMS']['job_ended'] = str(proc_pt_datetime_ended)


    # Write product config file for job.

    with open(product_config_ini_filename, 'w') as product_configfile:

        product_configfile.write("#" + "\n")
        product_configfile.write("# " + proc_utc_datetime + "\n")
        product_configfile.write("#" + "\n")
        product_configfile.write("# Machine-generated by " + swname + "\n")
        product_configfile.write("#" + "\n")
        product_configfile.write("\n")

        product_config.write(product_configfile)




    # Code-timing benchmark.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds after writing product config file at pipeline end =",
        end_time_benchmark - start_time_benchmark)
    start_time_benchmark = end_time_benchmark


    # Code-timing benchmark overall.

    end_time_benchmark = time.time()
    print("Elapsed time in seconds to run one instance of science pipeline =",
        end_time_benchmark - start_time_benchmark_at_start)


    # Termination.

    terminating_exitcode = 0
    if not run_sfft_was_successful:
        terminating_exitcode = 4

    print("terminating_exitcode =",terminating_exitcode)


    # AWS Batch job should be successful whenever terminating_exitcode < 64.

    aws_batch_job_exitcode = 0

    if (terminating_exitcode >= 64):
        aws_batch_job_exitcode = terminating_exitcode

    print("aws_batch_job_exitcode =",aws_batch_job_exitcode)

    exit(aws_batch_job_exitcode)








    exit(0)





