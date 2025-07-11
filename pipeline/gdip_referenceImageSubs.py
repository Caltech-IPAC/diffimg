import csv
import re
import boto3
from botocore.exceptions import ClientError
from astropy.io import fits
import numpy as np

import modules.utils.gdip_subs as util


# Subs used by the general difference-image pipeline related to reference images and catalogs.


#####################################################################################
# Generate reference-image catalog and upload it to S3 bucket.
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
    sextractor_refimage_dict["sextractor_PARAMETERS_NAME".lower()] = cfg_path + "/rapidSexParamsRefImage.inp"
    sextractor_refimage_dict["sextractor_FILTER_NAME".lower()] = cfg_path + "/cdf/rapidSexRefImageFilter.conv"
    sextractor_refimage_dict["sextractor_STARNNW_NAME".lower()] = cfg_path + "/rapidSexRefImageStarGalaxyClassifier.nnw"
    sextractor_refimage_dict["sextractor_CATALOG_NAME".lower()] = filename_refimage_catalog
    sextractor_cmd = util.build_sextractor_command_line_args(sextractor_refimage_dict)
    exitcode_from_sextractor = util.execute_command(sextractor_cmd)


    return
