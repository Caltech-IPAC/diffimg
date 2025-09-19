from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt

hdul = fits.open('bkg_subbed_science_image.fits')
data = hdul[0].data
wcs = WCS(hdul[0].header)
hdul.close()

hdul = fits.open('mosaic_6deg_j00056_1asec_cutout_resampled_gainmatched_convolved.fits')
data2 = hdul[0].data
wcs2 = WCS(hdul[0].header)
hdul.close()

hdul = fits.open('diffimage_masked.fits')
data3 = hdul[0].data
wcs3 = WCS(hdul[0].header)
hdul.close()

hdul = fits.open('naive_diffimage_masked.fits')
data4 = hdul[0].data
wcs4 = WCS(hdul[0].header)
hdul.close()

# Example:
center_pixel = (268, 1668)
cutout_size_pixels = (200, 200)
cutout = Cutout2D(data, position=center_pixel, size=cutout_size_pixels, wcs=wcs)
cutout2 = Cutout2D(data2, position=center_pixel, size=cutout_size_pixels, wcs=wcs)
cutout3 = Cutout2D(data3, position=center_pixel, size=cutout_size_pixels, wcs=wcs)
cutout4 = Cutout2D(data4, position=center_pixel, size=cutout_size_pixels, wcs=wcs)


'''
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 2, 1, projection=cutout.wcs)
ax.imshow(cutout.data, origin='lower', cmap='gray_r', vmin=-20, vmax=20)
ax.imshow(cutout2.data, origin='lower', cmap='gray_r', vmin=-20, vmax=20)


# Add WCS axes labels and grid
ax.coords.grid(True, color='white', ls='solid')
ax.coords[0].set_axislabel('Right Ascension')
ax.coords[1].set_axislabel('Declination')



plt.title('FITS Image Cutout')
plt.show()
'''





# Create a figure with 2 rows and 2 columns of subplots
fig, axs = plt.subplots(2, 2)

ct = 'rainbow'
vmin = -20
vmax = 300

axs[0, 0].imshow(cutout.data, origin='lower', cmap=ct, vmin=vmin, vmax=vmax)
axs[0, 0].set_title('Science Image')

axs[0, 1].imshow(cutout2.data, origin='lower', cmap=ct, vmin=vmin, vmax=vmax)
axs[0, 1].set_title('Reference Image')

axs[1, 0].imshow(cutout3.data, origin='lower', cmap=ct, vmin=vmin, vmax=vmax)
axs[1, 0].set_title('ZOGY Difference Image')

axs[1, 1].imshow(cutout4.data, origin='lower', cmap=ct, vmin=vmin, vmax=vmax)
axs[1, 1].set_title('Naive Difference Image')

# Adjust layout to prevent overlapping titles/labels
plt.tight_layout()

#plt.title('FITS Image Cutout')
plt.show()

