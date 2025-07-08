from astropy.wcs import WCS
from astropy.io import fits


input_fits_file = 'ADP.2022-07-27T14_56_30.297.fits'
hdr_num = 1           # Second HDU

# Load a FITS file
with fits.open(input_fits_file) as hdul:
    hdr = hdul[hdr_num].header

    radecsys = hdr['RADECSYS']
    print("radecsys =",radecsys)
    hdr.remove('RADECSYS', remove_all=True)
    hdr['RADESYSA'] = radecsys

    hdr.remove('PROJP1', remove_all=True)
    hdr.remove('PROJP3', remove_all=True)
    hdr.remove('PROJP5', remove_all=True)

    w = WCS(hdr) # Initialize WCS object from FITS header

print(w)

print("CTYPE = ",w.wcs.crpix)

naxis1 = hdr['NAXIS1']
naxis2 = hdr['NAXIS2']

print("naxis1,naxis2 =",naxis1,naxis2)

crpix1 = w.wcs.crpix[0]
crpix2 = w.wcs.crpix[1]


# Example of converting pixel coordinates to celestial coordinates
pixel_x, pixel_y = crpix1, crpix2
celestial_coords = w.pixel_to_world(pixel_x, pixel_y)
print(f"Pixel ({pixel_x}, {pixel_y}) corresponds to {celestial_coords.ra.deg:.12f} RA and {celestial_coords.dec.deg:.12f} Dec.")



pixel_x, pixel_y = 0,0
celestial_coords = w.pixel_to_world(pixel_x, pixel_y)
print(f"Pixel ({pixel_x}, {pixel_y}) corresponds to {celestial_coords.ra.deg:.12f} RA and {celestial_coords.dec.deg:.12f} Dec.")
