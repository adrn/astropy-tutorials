# coding: utf-8

"""  """

from __future__ import division, print_function

# Standard library
import os, sys

# Third-party
import numpy as np
from scipy.optimize import leastsq
import emcee

from astropy import wcs
from astropy.io import fits

hdulist = fits.open("L5g_0355+11_Cruz09.fits")
hdu = hdulist[0]
w = wcs.WCS(hdu.header)

flux = hdu.data
wavelength = w.all_pix2world(np.arange(len(flux))[:,np.newaxis], 0)
wavelength = wavelength.reshape(len(wavelength))

fig,axes = subplots(2,1,figsize=(12,10))
axes[0].plot(wavelength, flux, drawstyle="steps")
axes[0].set_xlim(wavelength.min(), wavelength.max())
axes[0].set_ylabel("Flux [{0}]".format(hdu.header["BUNIT"]), fontsize=20)
axes[0].axvline(6708, color="red", linestyle="--")
axes[0].xaxis.set_ticklabels([])

# Zoomed in panel
axes[1].plot(wavelength, flux, drawstyle="steps")
axes[1].set_xlim(6650, 6800)
axes[1].set_ylim(-0.001, 0.01)
axes[1].set_xlabel(r"Wavelength [$\AA$]", fontsize=20)
axes[1].set_ylabel("Flux [{0}]".format(hdu.header["BUNIT"]), fontsize=20)

fig.subplots_adjust(hspace=0.05)

idx = (wavelength > 6650) & (wavelength < 6750)
Li_flux = flux[idx]
Li_wvln = wavelength[idx]

def gaussian(x, mu, sigma):
    return 1. / np.sqrt(2*np.pi) / sigma * np.exp(-(x-mu)**2 / (2*sigma**2))
    
def model1(p, x):
    """ One Gaussian + straight line model """
    return p[0]*gaussian(x, p[1], p[2]) + p[3]*x + p[4]

def model2(p, x):
    """ Two Gaussians + straight line model """
    return p[0]*gaussian(x,p[1],p[2]) + p[3]*gaussian(x,p[4],p[5]) + p[6]*x + p[7]

def error_function1(p, x, y):
    return model1(p, x) - y

def error_function2(p, x, y):
    return model2(p, x) - y

p_opt, ier = leastsq(error_function1, [-0.01, 6700., 5., 0., 0.], args=(Li_wvln, Li_flux))
p_opt2, ier = leastsq(error_function2, [-0.01, 6700., 5., 0, 6710., 5., 0., 0.], args=(Li_wvln, Li_flux))

figure(figsize=(10,8))
plot(Li_wvln, Li_flux, drawstyle="steps", color="k")
plot(Li_wvln, model1(p_opt, Li_wvln), color="b")
plot(Li_wvln, model2(p_opt2, Li_wvln), color="r")
xlabel(r"Wavelength [$\AA$]", fontsize=20)
