import csv
import re
from astropy.io import fits
import numpy as np

import modules.utils.gdip_subs as util


# Subs used by the general difference-image pipeline related to reference images and catalogs.


#####################################################################################
# Generate reference-image catalog.
#####################################################################################

def generateReferenceImageCatalog(filename_refimage_image,
                                  filename_refimage_uncert,
                                  cfg_path,
                                  sextractor_refimage_dict,
                                  filename_refimage_catalog):


    # Compute SExtractor catalog for reference image.

    sextractor_refimage_dict["sextractor_detection_image".lower()] = "None"
    sextractor_refimage_dict["sextractor_input_image".lower()] = filename_refimage_image
    sextractor_refimage_dict["sextractor_WEIGHT_IMAGE".lower()] = filename_refimage_uncert
    sextractor_refimage_dict["sextractor_PARAMETERS_NAME".lower()] = cfg_path + "/srcExtractParamsRefImage.inp"
    sextractor_refimage_dict["sextractor_FILTER_NAME".lower()] = cfg_path + "/srcExtractRefImageFilter.conv"
    sextractor_refimage_dict["sextractor_STARNNW_NAME".lower()] = cfg_path + "/srcExtractRefImageStarGalaxyClassifier.nnw"
    sextractor_refimage_dict["sextractor_CATALOG_NAME".lower()] = filename_refimage_catalog
    sextractor_cmd = util.build_sextractor_command_line_args(sextractor_refimage_dict)
    exitcode_from_sextractor = util.execute_command(sextractor_cmd)


    return


#####################################################################################
# Compute reference-image coverage map (post facto with limited information).
#####################################################################################


def compute_coverage_map(input_filename,hdu_index,output_filename,nframes):

    hdul = fits.open(input_filename)
    hdr = hdul[hdu_index].header
    data = hdul[hdu_index].data

    hdul.close()

    hdr["BUNIT"] = "counts"

    np_data = np.array(data)

    cov_map = np.where(np.isnan(np_data), 0.0, nframes)

    hdu_cov = fits.PrimaryHDU(header=hdr,data=cov_map.astype(np.float32))
    hdu_list_cov = []
    hdu_list_cov.append(hdu_cov)
    hdu_cov = fits.HDUList(hdu_list_cov)
    hdu_cov.writeto(output_filename,overwrite=True,checksum=True)

    return
