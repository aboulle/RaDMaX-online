# coding=utf8
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from tools import  *
import numpy as np


# Debye-Waller Spline à N abscisses
def f_DW_spline3_smooth(alt, dwp,t):
	w_DW_free = dwp[:]
	w_DW = np.array([1.0,1.0,1.0])
	w_DW = np.append(w_DW,w_DW_free)
	N_abscisses = len(w_DW) - 3.
	z = alt * N_abscisses / t
	index = 0
	DW = np.ones(len(z))
	for i in z:
		DW[index] = cubicSpline(i,w_DW)
		index = index + 1
	return DW

def f_DW_spline3_abrupt(alt, dwp,t):
	w_DW = dwp[:]
	N_abscisses = len(w_DW) - 3.
	z = alt * N_abscisses / t
	index = 0
	DW = np.ones(len(z))
	for i in z:
		DW[index] = cubicSpline(i,w_DW)
		index = index + 1
	return DW

def f_DW_pv(alt, pv_p, t):
	height = 1 - pv_p[0]
	loc = pv_p[1] * t
	fwhm1 = pv_p[2] * t
	fwhm2 = pv_p[3] * t
	eta1 = pv_p[4]
	eta2 = pv_p[5]
	bkg = 1 - pv_p[6]
	DW = zeros(len(alt))
	DW[(alt<=loc)] = 1. - f_pVoigt(alt[alt<=loc], [height, loc, fwhm1, eta1])
	DW[(alt>loc)] = 1. - ( f_pVoigt(alt[alt>loc], [height-bkg, loc, fwhm2, eta2]) + bkg)
        #for n in arange(len(alt)):
            #if alt[n]<=loc:
                #DW[n] = 1. - f_pVoigt4nb(alt[n], array([height, loc, fwhm1, eta1],dtype=float64))
            #else:
                #DW[n] = 1 - (f_pVoigt4nb(alt[n], array([height-bkg, loc, fwhm2, eta2],dtype=float64)) + bkg)

	return DW

def control_dwp(alt,dwp,t,model):
	Nspline = len(dwp)
	Nslice = len(alt)-1

	if model == 'B-splines smooth':
		z_dwp =  np.arange(1, Nspline+1)* t / Nspline # generate depth (x axis) for the strain basis function
		z_interp = np.arange(1, (Nspline*Nslice)+1)* t / (Nspline*Nslice)
	if model == 'B-splines abrupt':
		z_dwp =  np.arange(0, Nspline)* t / (Nspline-1) # generate depth (x axis) for the strain basis function
		z_interp = np.arange(0, (Nspline*Nslice)* t / (Nspline*Nslice-1))
	dw = f_DW(z_interp, dwp, t, model)
	scaled_dwp = dw[ np.in1d(np.around(z_interp, decimals=1),  np.around(z_dwp,decimals=1))]

	return z_dwp, scaled_dwp

def f_DW(alt, dwp, t, model):
	if model == 'B-splines smooth':
		dw = f_DW_spline3_smooth(alt, dwp, t)
	if model == 'B-splines abrupt':
		dw = f_DW_spline3_abrupt(alt, dwp, t)
	return dw

def old2new_DW(alt, dwp, t, new_size,model):
     dwp_guess = np.ones(new_size)
     dw_old = f_DW(alt, dwp, t, model)
     def errfunc(dwp, alt, dw, t, model): return f_DW(alt, dwp, t, model) - dw_old
     dwp_new, success = leastsq(errfunc, dwp_guess, args=(alt, dw_old, t, model))
     return dwp_new

def interp_and_fit_dw(x, y, alt, old_dwp, model):
	# interpolate the dw curve obtained with the control points with a deg3 spline
	# fit a new a new dw curve to this to obtain new Bspline weights
	# old weights are used as a guess
	f = interp1d(x,y, kind = 'cubic', fill_value="extrapolate")
	int_curve = f(alt)
	t = alt.max()

	def errfunc(dwp, alt, int_curve):
		return f_DW(alt, dwp, t, model) - int_curve

	dwp_fit, success = leastsq(errfunc, old_dwp, args=(alt, int_curve))
	return dwp_fit

