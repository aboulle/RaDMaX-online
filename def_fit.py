# coding=utf8
import numpy as np
from scipy.optimize import least_squares
from def_XRD import *

def fit_curve(x, y, param, cst):
	Nw = cst["sdw_basis"]
	bounds_min = np.concatenate((np.full(Nw, cst["min_strain"]), np.full(Nw, cst["min_dw"])))
	bounds_max = np.concatenate((np.full(Nw, cst["max_strain"]), np.full(Nw, cst["max_dw"])))

	param[param<bounds_min] = bounds_min[param<bounds_min]
	param[param>bounds_max] = bounds_max[param>bounds_max]


	def residual(param, x, y, cst):
		y_cal = f_Refl(x, param, cst)[0]
		return (np.log10(y) - np.log10(y_cal))

	if cst["algo"] == 'Least squares':
		res_lsq = least_squares(residual, param, bounds=(bounds_min, bounds_max), args=(x, y,cst))
		param_fit = res_lsq.x
		return param_fit[0:Nw], param_fit[Nw: 2*Nw]

	if cst["algo"] == 'Least squares (no bounds)':
		res_lsq = least_squares(residual, param, method='lm', args=(x, y,cst))
		param_fit = res_lsq.x
		return param_fit[0:Nw], param_fit[Nw: 2*Nw]


	#t0 = time()

	#print ("Fitting... Please Wait")
	#par_fit, success = leastsq(residual, const.par, args = (const.Iobs, const.th))
	#print("Done")
	#print("Fitting time", "%5.4f" %(time() - t0), "secondes")
