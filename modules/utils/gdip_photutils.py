import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

import modules.utils.gdip_subs as util

from photutils.detection import DAOStarFinder
from photutils.psf import PSFPhotometry,ImagePSF,IterativePSFPhotometry


#-------------------------------------------------------------------
# Compute PSF-fit catalog with photutils.

def compute_diffimage_psf_catalog(n_clip_sigma,
                                  n_thresh_sigma,
                                  fwhm,
                                  fit_shape,
                                  aperture_radius,
                                  input_img_filename,
                                  input_unc_filename,
                                  input_psf_filename,
                                  output_psfcat_residual_filename):

    '''
    Method compute_diffimage_psf_catalog

    Inputs:
    n_clip_sigma            Number of sigmas for clipped data statistics.
    n_thresh_sigma          Number of sigmas above the background for source detection.
    fwhm                    Full-width half-maximum (FWHM) of Gaussian-kernel major axis, in units of pixels.
    fit_shape               Two-element tuple defining central star region, in pixels,
                            with highest flux signal-to-noise to be fit using the model PSF
    aperture_radius         Aperture radius to estimate the initial flux of each source.

    The threshold is the absolute image value above which to select sources, cast in terms of
    n_thresh_sigma times clipped background standard deviation, the same units as the input data image.

    Note that the input image must be background-subtracted prior to using the photometry classes.

    In general, fit_shape should be set to a small size (e.g., (5, 5)) that covers the central star region.

    Odd numbers for NAXIS1 and NAXIS2 is a required assumption for input PSF.
    '''


    print ("n_thresh_sigma =",n_thresh_sigma)
    print ("n_clip_sigma =",n_clip_sigma)
    print ("fwhm =",fwhm)
    print ("fit_shape =",fit_shape)
    print ("aperture_radius =",aperture_radius)
    print("input_img_filename =",input_img_filename)
    print("input_unc_filename =",input_unc_filename)
    print("input_psf_filename =",input_psf_filename)
    print("output_psfcat_residual_filename =",output_psfcat_residual_filename)

    hdu_index = 0
    saturation_level_image_rate = 999999
    stats_image = util.fits_data_statistics_with_clipping(input_img_filename,\
                                                          n_clip_sigma,\
                                                          hdu_index,\
                                                          saturation_level_image_rate)

    avg_image = stats_image["clippedavg"]
    std_image = stats_image["clippedstd"]
    cnt_image = stats_image["nkept"]
    noutliers_image = stats_image["noutliers"]
    gmed_image = stats_image["gmed"]
    datascale_image = stats_image["gsigma"]
    gmin_image = stats_image["gdatamin"]
    gmax_image = stats_image["gdatamax"]
    npixsat_image = stats_image["satcount"]
    npixnan_image = stats_image["nancount"]

    print("Image-data statistics: clippedavg, gmed, clippedstd = ",np.array((avg_image, gmed_image, std_image)))

    threshold = n_thresh_sigma * std_image
    print ("threshold =",threshold)

    hdul_image = fits.open(input_img_filename)
    hdr_image = hdul_image[0].header
    data_image = hdul_image[0].data

    hdul_image.close()

    hdul_uncert = fits.open(input_unc_filename)
    hdr_uncert = hdul_uncert[0].header
    data_uncert = hdul_uncert[0].data

    hdul_uncert.close()

    hdul_psf = fits.open(input_psf_filename)
    hdr_psf = hdul_psf[0].header
    data_psf = hdul_psf[0].data

    hdul_psf.close()

    naxis1 = hdr_psf["NAXIS1"]
    naxis2 = hdr_psf["NAXIS2"]

    print("naxis1,naxis2 =",naxis1,naxis2)


    # x_0 and y_0 is zero-based pixel index coordinates of PSF center.

    x_0 = (naxis1 - 1) / 2
    y_0 = (naxis2 - 1) / 2

    print("x_0,y_0 =", x_0,y_0)
    print("Type of x_0,y_0 =", type(x_0),type(y_0))

    data_psf_np = np.array(data_psf)
    data_psf_sum = np.sum(data_psf_np)

    print("data_psf_sum =", data_psf_sum)

    psf_model = ImagePSF(data_psf_np,flux=data_psf_sum,x_0=x_0,y_0=y_0)


    # Optionally plot the model PSF.

    if plot_flag:
        yy, xx = np.mgrid[:naxis1, :naxis2]
        data = psf_model(xx, yy)
        plt.imshow(data, origin='lower')
        plt.show()


    # Initialize the DAOStarFinder and PSFPhotometry class instances.

    finder = DAOStarFinder(threshold, fwhm)

    finder_attributes = finder.__dict__.keys()

    print("finder_attributes =",finder_attributes)


    # iterative_flag must be false for pipeline, but can be True temporarily for experiments.
    # Default fitter_maxiters=100 iterations seems to add artifacts to the residual image.

    iterative_flag = False

    if not iterative_flag:
        try:
            psfphot = PSFPhotometry(psf_model=psf_model,fit_shape=fit_shape,finder=finder,aperture_radius=aperture_radius)
            psfcat_flag = True
        except:
            print("*** Warning: Could not make psf-fit catalog (perhaps no sources were detected); continuing...")
            psfcat_flag = False
            psfphot = None
    else:
        psfphot = IterativePSFPhotometry(psf_model=psf_model,fit_shape=fit_shape,finder=finder,aperture_radius=aperture_radius)

    psfphot_attributes = psfphot.__dict__.keys()

    print("psfphot_attributes =",psfphot_attributes)


    # Call PSFPhotometry class instance on the data array to do the PSF-fitting to the image data.

    try:
        phot = psfphot(data=data_image,error=data_uncert,filter_non_finite=True)
    except:
        print("*** Warning: Exception thrown calling PSFPhotometry class instance on the data array to do the PSF-fitting to the image data; continuing...")
        psfcat_flag = False
        phot = None


    # Make residual image.

    if psfcat_flag:
        try:
            resid = psfphot.make_residual_image(data_image)
            new_hdu = fits.PrimaryHDU(data=resid.astype(np.float32))
            new_hdu.writeto(output_psfcat_residual_filename,overwrite=True,checksum=True)
        except:
            print("*** Warning: Could not make residual image (perhaps no sources were detected); continuing...")


    # Return photometry and finder catalogs.

    return psfcat_flag,phot,psfphot
