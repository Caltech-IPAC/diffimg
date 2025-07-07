import pyvo as vo
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.utils.data import download_file
from astropy.io import fits
import re
import os


# Specify sky position of interest and desired filter.

ra = 314.30417
dec = 77.595559
pos = SkyCoord(ra=ra, dec=dec, unit='deg')
2mass_filter = 'H'

# Query for 2MASS images that overlap sky position within 1.0 arcseconds.

twomass_service = vo.dal.SIAService("https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?type=at&ds=asky&")
im_table = twomass_service.search(pos=pos, size=1.0*u.arcsec)
print(im_table.to_table())


# Columns in table.
# name,download,center_ra,center_dec,naxes,naxis,scale,format,crpix,crval,crota2,band,bref,bhi,blo,pers_art,glint_art,type,dataset,pixflags,id,scntr,date,hem,scan,image,ut_date,coadd_key,seesh,magzp,msnr10,bin


# Get the first returned 2MASS image for the filter of interest.

for i in range(len(im_table)):
    if im_table[i]['band'] == 2mass_filter:
        break
print(im_table[i].getdataurl())


# Download the 2MASS image and decompress it.
#
# curl --output hi0600256.fits.gz "https://irsa.ipac.caltech.edu:443/cgi-bin/2MASS/IM/nph-im?ds=asky&atdir=/ti08&dh=000616n&scan=060&name=hi0600256.fits"
# gunzip hi0600256.fits.gz

download_url = im_table[i].getdataurl()

filename_match = re.match(r".+?name\=(.+?.fits)", download_url)

try:
    fits_file = filename_match.group(1)
    print("fits_file =",fits_file)

    gz_fits_file = fits_file + ".gz"

    curl_cmd = "curl --output " + gz_fits_file + " \"" + download_url + "\""
    print("curl_cmd =",curl_cmd)

    return_code = os.system(curl_cmd)
    print(f"Command exited with code: {return_code}")

    gunzip_cmd = "gunzip " + gz_fits_file
    print("gunzip_cmd =",gunzip_cmd)

    return_code = os.system(gunzip_cmd)
    print(f"Command exited with code: {return_code}")

except:
    print("*** Error: No FITS filename match found; quitting...")
    exit(64)

exit(0)



