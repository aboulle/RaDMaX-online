# coding=utf8
import numpy as np
from scipy.optimize import least_squares, dual_annealing
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

	if cst["algo"] == 'Test':
		res_lsq = dual_annealing(residual, x0=param, bounds=np.c_[bounds_min, bounds_max], restart_temp_ratio=0.0, args=(x, y,cst))
		param_fit = res_lsq.x
		return param_fit[0:Nw], param_fit[Nw: 2*Nw]
	
	if cst["algo"] == 'Least squares':
		res_lsq = least_squares(residual, param, bounds=(bounds_min, bounds_max), args=(x, y,cst))
		param_fit = res_lsq.x
		return param_fit[0:Nw], param_fit[Nw: 2*Nw]

	if cst["algo"] == 'Least squares (no bounds)':
		res_lsq = least_squares(residual, param, method='lm', args=(x, y,cst))
		param_fit = res_lsq.x
		return param_fit[0:Nw], param_fit[Nw: 2*Nw]

def SMB_fit(residual, x, y, cst, p0, bounds=(-np.inf, np.inf), max_nfev=2, mbmult = 1, epochs = 100):
	
	for epoch in range(epochs):
		print(epoch)
		m = len(x)
		minibatch_size = mbmult * len(p0)
		
		shuffled_indices = np.random.permutation(m)
		x_shuffled = x[shuffled_indices]
		y_shuffled = y[shuffled_indices]
		
		for i in range(0, m, minibatch_size):
			xi = x[i:i+minibatch_size]
			yi = y[i:i+minibatch_size]
			
			cst["resol"] = f_resol(xi, p0, "Pseudo-Voigt")
			
			mb_lsq = least_squares(residual, p0, bounds=bounds, max_nfev=max_nfev, args=(xi, yi, cst))
			p0 = mb_lsq.x
	
	return p0

	

	#t0 = time()

	#print ("Fitting... Please Wait")
	#par_fit, success = leastsq(residual, const.par, args = (const.Iobs, const.th))
	#print("Done")
	#print("Fitting time", "%5.4f" %(time() - t0), "secondes")
