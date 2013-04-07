

from __future__ import division, print_function

# Standard library
import os, sys

# Third-party
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import leastsq
import emcee

#Astropy
from astropy import wcs
from astropy.io import fits

data_file = 'L5g_0355+11_Cruz09.fits'

hdulist = fits.open(data_file) #read in all HDUs/extensions
#in this case, only one HDU
hdu0 = hdulist[0]
wcs0 = wcs.WCS(hdu0.header)

#in this case, the only data in the FITS file is the flux values
flux0 = hdu0.data
n_pix0=len(flux0)

# create numpy array of integers from 0 to number of pixels as a row. then reshaped to be column. placeholder for wavelength. 
# IDL equivalent = INTARR(n_elements(flux))
pixel_grid = np.arange(0,n_pix0).reshape(n_pix0,1)

#make wavelength array with astropy.wcs magic via all_pix2world method
wavelength0 = wcs0.all_pix2world(pixel_grid,0)

#reshape to row
wavelength0 = wavelength0.reshape(n_pix0)

##Now we have wavelength and flux for 0th HDU/extension.


# measure EW of Lithium line.
idx = (wavelength0 > 6690) & (wavelength0 < 6725)
Li_flux = flux0[idx]
Li_wvln = wavelength0[idx]

# create gaussian function (this should be built-in to astropy)
def gaussian(x, mu, sigma):
    return 1. / np.sqrt(2*np.pi) / sigma * np.exp(-(x-mu)**2 / (2*sigma**2))

def model(p, x):
    """ Gaussian + constant line """
    return p[0]*gaussian(x, p[1], p[2]) + p[3]

def error_function(p, x, y):
    """ Generic error function given a model and data, assuming 
        no uncertainties 
    """
    return model(p, x) - y

leastsq_parameters, ier = leastsq(error_function, 
                            [-0.01, 6705., 5., 0.], 
                            args=(Li_wvln, Li_flux))

fig,axes = plt.subplots(2,1,figsize=(12,10))
axes[0].plot(wavelength0, flux0, drawstyle="steps")
axes[0].set_xlim(wavelength0.min(), wavelength0.max())
axes[0].set_ylabel("Flux [{0}]".format(hdu0.header["BUNIT"]), fontsize=20)
axes[0].axvline(6708, color="red", linestyle="--")
axes[0].xaxis.set_ticklabels([])

# Zoomed in panel
axes[1].plot(wavelength0, flux0, drawstyle="steps")
axes[1].set_xlim(6650, 6800)
axes[1].set_ylim(-0.001, 0.01)
axes[1].set_xlabel(r"Wavelength [$\AA$]", fontsize=20)
axes[1].set_ylabel("Flux [{0}]".format(hdu0.header["BUNIT"]), fontsize=20)

fig.subplots_adjust(hspace=0.05)

fig,axes = plt.subplots(2,1,figsize=(10,12))
fig.suptitle("Using scipy.optimize.leastsq")

# plot the individual model components
axes[0].plot(Li_wvln, leastsq_parameters[0]*gaussian(Li_wvln, *leastsq_parameters[1:3]))
axes[0].plot(Li_wvln, np.ones_like(Li_wvln)*leastsq_parameters[3])
axes[0].set_xticklabels([])

# plot the sum of the model components over the data
axes[1].plot(Li_wvln, Li_flux, drawstyle="steps", color="k")
print("Line center: {0:.2f}".format(leastsq_parameters[1]))
axes[1].plot(Li_wvln, model(leastsq_parameters, Li_wvln), color="b", linewidth=2)
axes[1].set_xlabel(r"Wavelength [$\AA$]", fontsize=20)
fig.subplots_adjust(hspace=0.08)

plt.show()