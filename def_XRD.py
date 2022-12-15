"""
Copyright (C) 2022 CNRS
Author(s) Alexandre BOULLE
…
This software is governed by the CeCILL Free Software License Agreement v2.1
You can  use, modify and/ or redistribute the software under the terms of the CECILL-2.1 at the following URL https://spdx.org/licenses/CECILL-2.1.html
The fact that you are presently reading this means that you have had knowledge of the CECILL-2.1 and that you accept its terms.
"""
##coding: Utf-8
from def_strain import *
from def_DW import *
#from tools import *
from time import *
import numpy as np

try:
    import TakagiTaupin as TT
    fast = True
except:
    pass

# compute list of constants used in the XRD computation
def compute_cst(w_list, th):
	re = 2.818*10**-5
	irf = w_list[0].value
	if irf == 'Split pseudo-Voigt':
		width = w_list[1].value
		try:
			a,b = width.split(',')
		except:
			a,b = float(w_list[1].value),float(w_list[1].value)
		width = (float(a)*np.pi/180,float(b)*np.pi/180)

		eta = w_list[2].value
		try:
			a,b = eta.split(',')
		except:
			a,b = float(w_list[2].value),float(w_list[2].value)
		eta = (float(a),float(b))

	else:
		try:
			width = float(w_list[1].value)*np.pi/180
		except:
			width = float((w_list[1].value).split(',')[0])*np.pi/180
		try:
			eta = float(w_list[2].value)
		except:
			eta = float((w_list[2].value).split(',')[0])



	wl =  float(w_list[3].value)
	offset = float(w_list[4].value)
	bkg = float(w_list[5].value)

	mater = xu.materials.Crystal.fromCIF(os.path.join("structures", w_list[6].value))
	a = mater.lattice.a
	b = mater.lattice.b
	c = mater.lattice.c
	alpha = mater.lattice.alpha
	beta = mater.lattice.beta
	gamma = mater.lattice.gamma
	Vol = mater.lattice.UnitCellVolume()
	h, k, l = int(w_list[7].value), int(w_list[8].value), int(w_list[9].value)
	phi = 0. #TODO implement asymmetric scattering
	d = mater.planeDistance(h,k,l)
	FH = mater.StructureFactor(mater.Q(h,k,l),12398/wl)
	FmH = mater.StructureFactor(mater.Q(-h,-k,-l),12398/wl)
	F0 = mater.StructureFactor(mater.Q(0,0,0),12398/wl)

	G = re * wl * wl / (np.pi * Vol) # Gamma : Structure factor -> Polarizability
	thB_S  = np.arcsin(wl / (2*d))
	g0 = np.sin(thB_S - phi) # gamma 0
	gH = -np.sin(thB_S + phi) # gamma H
	b_S = g0 / gH

	t = float(w_list[13].value)*10
	N = float(w_list[14].value) # number of slices
	t_l = t/N # slice thickness
	z = np.arange(N+1) * t_l # distance from interface

	resol = f_resol(th, [1, (th.min()+th.max())/2 , width, eta], irf)
	sample = w_list[10].value
	sdw_model = w_list[11].value
	sdw_basis = int(w_list[12].value) # number of splines

	algo = w_list[20].value
	min_strain = float(w_list[21].value)
	max_strain = float(w_list[22].value)
	min_dw = float(w_list[23].value)
	max_dw = float(w_list[24].value)
	gsa_temp = float(w_list[25].value)
	gsa_cycles = float(w_list[26].value)
	gsa_Tsteps = float(w_list[27].value)

	cst = {"re":re, "wl":wl, "bkg":bkg, "offset":offset, "width":width, "eta":eta,
			"a": a, "b":b, "c":c, "alpha":alpha, "beta":beta, "gamma":gamma,
			"Vol":Vol, "phi":phi, "FH":FH, "FmH":FmH, "F0":F0,
			"G":G, "thB_S":thB_S, "g0":g0, "gH":gH, "b_S":b_S,
			"t":t, "N":N, "t_l":t_l, "z":z, "resol":resol,
			"sample":sample, "sdw_model":sdw_model, "sdw_basis":sdw_basis,
			"algo":algo, "min_strain":min_strain, "max_strain":max_strain,
			 "min_dw":min_dw, "max_dw":max_dw, "gsa_temp":gsa_temp,
			 "gsa_cycles":gsa_cycles, "gsa_Tsteps":gsa_Tsteps}

	if (sample == 'Single crystal') or (sample =='Thin film'):
		return cst

	elif (sample == 'Thick film'):
		cst["t_film"] = int(w_list[15].value)*10
		return cst

	else:
		cst["t_film"] = int(w_list[15].value)*10
		sub = xu.materials.Crystal.fromCIF(os.path.join("structures", w_list[16].value))
		cst["a_sub"] = sub.lattice.a
		cst["b_sub"] = sub.lattice.b
		cst["c_sub"] = sub.lattice.c
		cst["alpha_sub"] = sub.lattice.alpha
		cst["beta_sub"] = sub.lattice.beta
		cst["gamma_sub"] = sub.lattice.gamma
		cst["Vol_sub"] = sub.lattice.UnitCellVolume()
		h_sub, k_sub, l_sub = int(w_list[17].value), int(w_list[18].value), int(w_list[19].value)
#		cst["h_sub"], cst["k_sub"], cst["l_sub"] =  h_sub, k_sub, l_sub
		cst["phi_sub"] = 0. #TODO implement asymmetric scattering
		d_sub = sub.planeDistance(h_sub,k_sub,l_sub)
		cst["FH_sub"] = sub.StructureFactor(sub.Q(h_sub,k_sub,l_sub),12398/wl)
		cst["FmH_sub"] = sub.StructureFactor(sub.Q(-h_sub,-k_sub,-l_sub),12398/wl)
		cst["F0_sub"] =sub.StructureFactor(sub.Q(0,0,0),12398/wl)

		cst["G_sub"] = re * wl * wl / (np.pi * cst["Vol_sub"]) # Gamma : Structure factor -> Polarizability
		thB_sub  = np.arcsin(wl / (2*d_sub))
		cst["thB_sub"]  = thB_sub
		cst["g0_sub"] = np.sin(thB_sub - phi) # gamma 0
		cst["gH_sub"] = -np.sin(thB_sub + phi) # gamma H
		cst["b_sub"] = cst["g0_sub"] / cst["gH_sub"]
		return cst

def f_Refl_Default(th, param, cst):
	offset = cst["offset"]*np.pi/360
	G = cst["G"]
	thB_S = cst["thB_S"]
	wl = cst["wl"]
	t = cst["t"]
	N = cst["N"]
	resol = cst["resol"]
	b_S = cst["b_S"]
	phi = cst["phi"]
	t_l = cst["t_l"]
	z = cst["z"]
	FH = cst["FH"]
	FmH = cst["FmH"]
	F0 = cst["F0"]
	Nspline = cst["sdw_basis"]
	model = cst["sdw_model"]
	bkg = cst["bkg"]

	th = th + offset
	strain = f_strain(z, param[:Nspline:], t, model)
	DW = f_DW(z, param[Nspline:2*Nspline:], t, model)
	thB = thB_S - strain * np.tan(thB_S)

	eta = (-b_S*(th-thB_S)*np.sin(2*thB_S) - 0.5*G*F0*(1-b_S)) / ((abs(b_S)**0.5)*G*(FH*FmH)**0.5 )
	res = (eta - np.sign(eta.real)*((eta*eta - 1)**0.5))

	n = 1
	while (n<=N):
		g0 = np.sin(thB[n] - phi)
		gH = -np.sin(thB[n] + phi)
		b = g0 / gH
		T = np.pi * G * ((FH*FmH)**0.5) * t_l * DW[n]/ (wl * (abs(g0*gH)**0.5) )
		eta = (-b*(th-thB[n])*np.sin(2*thB_S) - 0.5*G*F0*(1-b)) / ((abs(b)**0.5)*G*DW[n]*(FH*FmH)**0.5)
		sqrt_eta2 = (eta*eta-1)**0.5

		S1 = (res - eta + sqrt_eta2)*np.exp(-1j*T*sqrt_eta2)
		S2 = (res - eta - sqrt_eta2)*np.exp(1j*T*sqrt_eta2)

		res = (eta + sqrt_eta2*((S1+S2)/(S1-S2)))
		n += 1

	ical = np.convolve(abs(res)**2, resol, mode='same')
	return (ical/ical.max())+bkg

def f_Refl_Thin_Film(th, param, cst):
	offset = cst["offset"]*np.pi/360
	G = cst["G"]
	thB_S = cst["thB_S"]
	wl = cst["wl"]
	t = cst["t"]
	N = cst["N"]
	resol = cst["resol"]
	b_S = cst["b_S"]
	phi = cst["phi"]
	t_l = cst["t_l"]
	z = cst["z"]
	FH = cst["FH"]
	FmH = cst["FmH"]
	F0 = cst["F0"]
	Nspline = cst["sdw_basis"]
	model = cst["sdw_model"]
	bkg = cst["bkg"]

	th = th + offset
	strain = f_strain(z, param[:Nspline:], t, model)
	DW = f_DW(z, param[Nspline:2*Nspline:], t, model)
	thB = thB_S - strain * np.tan(thB_S) ## angle de Bragg dans chaque lamelle

	eta = 0
	res = 0

	n = 1
	while (n<=N):
		g0 = np.sin(thB[n] - phi)
		gH = -np.sin(thB[n] + phi)
		b = g0 / gH
		T = np.pi * G * ((FH*FmH)**0.5) * t_l * DW[n]/ (wl * (abs(g0*gH)**0.5) )
		eta = (-b*(th-thB[n])*np.sin(2*thB_S) - 0.5*G*F0*(1-b)) / ((abs(b)**0.5)*G*DW[n]*(FH*FmH)**0.5)
		sqrt_eta2 = (eta*eta-1)**0.5

		S1 = (res - eta + sqrt_eta2)*np.exp(-1j*T*sqrt_eta2)
		S2 = (res - eta - sqrt_eta2)*np.exp(1j*T*sqrt_eta2)

		res = (eta + sqrt_eta2*((S1+S2)/(S1-S2)))
		n += 1

	ical = np.convolve(abs(res)**2, resol, mode='same')
	return (ical/ical.max())+bkg

def f_Refl_Thick_Film(th, param,cst):
	offset = cst["offset"]*np.pi/360
	G = cst["G"]
	thB_S = cst["thB_S"]
	wl = cst["wl"]
	#h, k, l =cst["h"], cst["k"], cst["l"]
	t = cst["t"]
	N = cst["N"]
	resol = cst["resol"]
	b_S = cst["b_S"]
	phi = cst["phi"]
	t_l = cst["t_l"]
	t_film = cst["t_film"]
	z = cst["z"]
	FH = cst["FH"]
	FmH = cst["FmH"]
	F0 = cst["F0"]
	Nspline = cst["sdw_basis"]
	model = cst["sdw_model"]
	bkg = cst["bkg"]

	delta_t = t_film - t

	th = th + offset
	strain = f_strain(z, param[:Nspline:], t, model)
	DW = f_DW(z, param[Nspline:2*Nspline:], t, model)
	thB = thB_S - strain * np.tan(thB_S)

	eta = 0
	res = 0

	g0 = np.sin(thB[0] - phi)
	gH = -np.sin(thB[0] + phi)
	b = g0 / gH
	T = np.pi * G * ((FH*FmH)**0.5) * delta_t / (wl * (abs(g0*gH)**0.5) )
	eta = (-b*(th-thB[0])*np.sin(2*thB_S) - 0.5*G*F0*(1-b)) / ((abs(b)**0.5)*G*(FH*FmH)**0.5)
	sqrt_eta2 = (eta*eta-1)**0.5

	S1 = (res - eta + sqrt_eta2)*np.exp(-1j*T*sqrt_eta2)
	S2 = (res - eta - sqrt_eta2)*np.exp(1j*T*sqrt_eta2)

	res = (eta + sqrt_eta2*((S1+S2)/(S1-S2)))

	n = 1
	while (n<=N):
		g0 = np.sin(thB[n] - phi)
		gH = -np.sin(thB[n] + phi)
		b = g0 / gH
		T = np.pi * G * ((FH*FmH)**0.5) * t_l * DW[n]/ (wl * (abs(g0*gH)**0.5) )
		eta = (-b*(th-thB[n])*np.sin(2*thB_S) - 0.5*G*F0*(1-b)) / ((abs(b)**0.5)*G*DW[n]*(FH*FmH)**0.5)
		sqrt_eta2 = (eta*eta-1)**0.5

		S1 = (res - eta + sqrt_eta2)*np.exp(-1j*T*sqrt_eta2)
		S2 = (res - eta - sqrt_eta2)*np.exp(1j*T*sqrt_eta2)

		res = (eta + sqrt_eta2*((S1+S2)/(S1-S2)))
		n += 1

	ical = np.convolve(abs(res)**2, resol, mode='same')
	return (ical/ical.max())+bkg

def f_Refl_Thick_Film_and_Substrate(th, param,cst):
	offset = cst["offset"]*np.pi/360
	G = cst["G"]
	thB_S = cst["thB_S"]
	wl = cst["wl"]
	t = cst["t"]
	N = cst["N"]
	resol = cst["resol"]
	b_S = cst["b_S"]
	phi = cst["phi"]
	t_l = cst["t_l"]
	t_film = cst["t_film"]
	z = cst["z"]
	FH = cst["FH"]
	FmH = cst["FmH"]
	F0 = cst["F0"]
	Nspline = cst["sdw_basis"]
	model = cst["sdw_model"]
	bkg = cst["bkg"]
	thB_sub = cst["thB_sub"]
	b_sub = cst["b_sub"]
	F0_sub = cst["F0_sub"]
	FH_sub = cst["FH_sub"]
	FmH_sub = cst["FmH_sub"]

	delta_t = t_film - t

	th = th + offset
	strain = f_strain(z, param[:Nspline:], t, model)
	DW = f_DW(z, param[Nspline:2*Nspline:], t, model)
	thB = thB_S - strain * np.tan(thB_S)

	eta = (-b_S*(th-thB_sub)*np.sin(2*thB_sub) - 0.5*G*F0_sub*(1-b_sub)) / ((abs(b_sub)**0.5)*G*(FH_sub*FmH_sub)**0.5 )
	res = (eta - np.sign(eta.real)*((eta*eta - 1)**0.5))

	g0 = np.sin(thB[0] - phi)
	gH = -np.sin(thB[0] + phi)
	b = g0 / gH
	T = np.pi * G * ((FH*FmH)**0.5) * delta_t / (wl * (abs(g0*gH)**0.5) )
	eta = (-b*(th-thB[0])*np.sin(2*thB_S) - 0.5*G*F0*(1-b)) / ((abs(b)**0.5)*G*(FH*FmH)**0.5)
	sqrt_eta2 = (eta*eta-1)**0.5

	S1 = (res - eta + sqrt_eta2)*np.exp(-1j*T*sqrt_eta2)
	S2 = (res - eta - sqrt_eta2)*np.exp(1j*T*sqrt_eta2)

	res = (eta + sqrt_eta2*((S1+S2)/(S1-S2)))

	n = 1
	while (n<=N):
		g0 = np.sin(thB[n] - phi)
		gH = -np.sin(thB[n] + phi)
		b = g0 / gH
		T = np.pi * G * ((FH*FmH)**0.5) * t_l * DW[n]/ (wl * (abs(g0*gH)**0.5) )
		eta = (-b*(th-thB[n])*np.sin(2*thB_S) - 0.5*G*F0*(1-b)) / ((abs(b)**0.5)*G*DW[n]*(FH*FmH)**0.5)
		sqrt_eta2 = (eta*eta-1)**0.5

		S1 = (res - eta + sqrt_eta2)*np.exp(-1j*T*sqrt_eta2)
		S2 = (res - eta - sqrt_eta2)*np.exp(1j*T*sqrt_eta2)

		res = (eta + sqrt_eta2*((S1+S2)/(S1-S2)))
		n += 1

	ical = np.convolve(abs(res)**2, resol, mode='same')
	return (ical/ical.max())+bkg



def f_Refl_choice(th, param, cst):
	if cst["sample"] == 'Single crystal':
		return f_Refl_Default(th, param, cst)
	if cst["sample"] == 'Thin film':
		return f_Refl_Thin_Film(th, param, cst)
	if cst["sample"] == 'Thick film':
		return f_Refl_Thick_Film(th, param, cst)
	if cst["sample"] == "Thick film + substrate":
		return f_Refl_Thick_Film_and_Substrate(th, param,cst)

def f_Refl(th, param, cst):
	t0 = time()
	ical = f_Refl_choice(th, param, cst)
	t = time() - t0
	return ical, t
