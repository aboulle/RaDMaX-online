# coding=utf8
#import csv
#from scipy import *
import os
import numpy as np
import xrayutilities as xu
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
#import math


# Find the coordinate of a value within array
def find(x, Z):
	i = where(x == Z)
	return i[0]

# sign function
def signe(x):
	return x.real / abs(x.real)

def f_pVoigt(x, param):
	max = param[0]
	pos = param[1]
	FWHM = param[2]
	eta = param[3]
	gauss = max * np.exp(-np.log(2.) * ((x-pos)/(0.5*FWHM))**2)
	lorentz = max / (1. + ((x - pos)/(0.5*FWHM))**2)
	return eta*lorentz + (1-eta)*gauss

def f_Gauss(x, param):
	max = param[0]
	pos = param[1]
	FWHM = param[2]
	gauss = max * np.exp(-np.log(2.) * ((x-pos)/(0.5*FWHM))**2)
	return gauss

def f_Lorentz(x, param):
	max = param[0]
	pos = param[1]
	FWHM = param[2]
	lorentz = max / (1. + ((x - pos)/(0.5*FWHM))**2)
	return lorentz

def f_gbell(x, param):
    max = param[0]
    pos = param[1]
    FWHM = param[2]
    b = param[3]
    return max / (1. + (abs((x - pos)/(0.5*FWHM)))**(2*b))

def f_splitpV(x, param):
    max_ = param[0]
    pos = param[1]
    FWHMl = param[2][0]
    FWHMr = param[2][1]
    etal = param[3][0]
    etar = param[3][1]
    pv = f_pVoigt(x, [max_, pos, FWHMl, etal])
    pvr = f_pVoigt(x, [max_, pos, FWHMr, etar])
    pv[x > pos] = pvr[x > pos]
    return pv

def f_resol(x, param, irf):
	if irf == "Pseudo-Voigt":
		return f_pVoigt(x, param)
	if irf == "Gaussian":
		return f_Gauss(x, param)
	if irf == "Lorentzian":
		return f_Lorentz(x, param)
	if irf == "Generalized bell":
		return f_gbell(x,param)
	if irf == "Split pseudo-Voigt":
		return f_splitpV(x,param)


# B-spline basis function
def bSpline3(z):
    if z <= 0 :
        return 0
    elif z <= 1:
        return (z**3)/6
    elif z <= 2:
        return (2./3) + z * (-2 + z*(2 - z/2))
    elif z <= 3:
        return (-22./3) + z * (10 + z*(-4 + z/2))
    elif z <= 4:
        return (32./3) + z * (-8 + z*(2 - z/6))
    else :
        return 0

# Cubic spline
def cubicSpline(z,w):
    somme = 0
    index = 0
    for poids in w:
        somme = somme + poids * bSpline3(z-index+3)
        index = index + 1
    return somme

# Compute root mean-squared error of log I
def rmse(yo, yc):
    mse = ((np.log10(yo)-np.log10(yc))**2).sum()/len(yo)
    return np.sqrt(mse)

def auto_strain(tth, iexp, cst):
    wl = cst["wl"]

    iexp[iexp==0]=1
    iexp /= iexp.max()
    Q = 4 * np.pi * np.sin(tth*np.pi/360)/wl 
    eps = 100 * (Q - Q[iexp==iexp.max()]) / Q

    ifilt = gaussian_filter1d(iexp, 2)
    peaks, _ = find_peaks(ifilt,
                               #height = 0,
                               threshold=1e-6,
                               distance = 1,
                               prominence = 1e-5)
    return -eps[peaks[0]]

