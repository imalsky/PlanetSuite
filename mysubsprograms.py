#!/usr/bin/env python

import math
import numpy as np
import os
import shutil
import scipy
from scipy import loadtxt, optimize
import sys
import time
import random
import pandas as pd
from scipy.interpolate import interp1d
from scipy import interpolate


msun = 1.9892e33
rsun = 6.9598e10
rearth = 6.371008e8
mjup = 1.8986e30
rjup = 6.9911e9
mearth = 5.97e27
sigma=5.67e-5
au = 1.496e13


def calculate_rho(mp, enFrac):
	observed_Mcore, observed_Rcore = loadtxt('coreMRcomp2_v40_all.txt', unpack=True, skiprows =11, usecols=[0,1])
	core_radius_function = interp1d(observed_Mcore, observed_Rcore, fill_value="extrapolate")

	#In units of Rearth and Mearth
	planet_core_mass = (mp * (1.0 - enFrac))

	planet_core_radius = core_radius_function(planet_core_mass)


	core_mass_cgs = planet_core_mass * mearth
	core_radius_cgs = planet_core_radius * rearth

	core_volume = (4./3.) * (3.14159) * (core_radius_cgs ** 3.0)
	rhocore = core_mass_cgs / float(core_volume)
	
	return (rhocore, planet_core_radius)


def calculate_column_depth(Teq, profile, Teff_star):
    e = 2.7182818284

    if Teff_star < 3500:
        T,P, k_r, k_p = loadtxt('OpacityTableSolarMetal.txt', unpack=True, skiprows =38, usecols=[0,1,5,6])
    elif 3500 <= Teff_star < 4500:
        T,P, k_r, k_p = loadtxt('OpacityTableSolarMetal.txt', unpack=True, skiprows =38, usecols=[0,1,7,8])
    elif 4500 <= Teff_star < 5500:
        T,P, k_r, k_p = loadtxt('OpacityTableSolarMetal.txt', unpack=True, skiprows =38, usecols=[0,1,9,10])
    elif 5500 <= Teff_star < 6500:
        T,P, k_r, k_p = loadtxt('OpacityTableSolarMetal.txt', unpack=True, skiprows =38, usecols=[0,1,11,12])
    else:
        T,P, k_r, k_p = loadtxt('OpacityTableSolarMetal.txt', unpack=True, skiprows =38, usecols=[0,1,13,14])

    Opacity_function = interpolate.interp2d(T, P, k_p)
    zone, mass, temperature, radius, pressure = loadtxt(profile, unpack=True, skiprows =6, usecols=[0,1,2,3,6])   

    switch_zone = []
    for i in range(len(zone)):
        mass_column_depth = ((mass[0] - mass[i]) * msun) / (4 * 3.14159 * (radius[i] ** 2))
        opacity_column_depth = (2 / (Opacity_function(Teq, pressure[i])))[0]
        switch_zone.append((opacity_column_depth - mass_column_depth, zone[i], mass_column_depth))

    column_depth = abs(switch_zone[0][0])
    
    if switch_zone[0][0] > 0:
        for i in range(len(switch_zone)):
            if switch_zone[i][0] < 0:
                column_depth = switch_zone[i - 1][2]
                break
    else:
        for i in range(len(switch_zone)):
            if switch_zone[i][0] > 0:
                column_depth = switch_zone[i - 1][2]
                break

    return column_depth



def run_pre_reduce(inlist_pre_reduce, initial_mod, pre_reduce_mod, mp):
	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_pre_reduce', 'r')
	g = f.read()
	f.close()

	g = g.replace("<<loadfile>>",'"' + initial_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + pre_reduce_mod + '"')
	g = g.replace("<<mp>>",str(5.0 * mearth / msun))

	h = open(inlist_pre_reduce, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_pre_reduce, "inlist")

	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time


def run_pre_core(inlist_pre_core, pre_reduce_mod, pre_core_mod, enFrac,core_mass,rho):
	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_pre_core', 'r')
	g = f.read()
	f.close()

	g = g.replace("<<loadfile>>",'"' + pre_reduce_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + pre_core_mod + '"')
	g = g.replace("<<core_mass>>", str(0.1 * core_mass * mearth / msun))
	g = g.replace("<<rho>>", str(rho))
	
	h = open(inlist_pre_core, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_pre_core, "inlist")

	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time



def run_comp(inlist_comp, pre_core_mod, comp_mod, z, y):
	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_comp', 'r')
	g = f.read()
	f.close()

	g = g.replace("<<pre_core_mod>>",'"' + pre_core_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + comp_mod + '"')
	g = g.replace("<<y>>",str(y))
	g = g.replace("<<z>>",str(z))

	h = open(inlist_comp, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_comp, "inlist")

	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time


def run_corel(inlist_corel, comp_mod, corel_mod):
	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_corel', 'r')
	g = f.read()
	f.close()

	g = g.replace("<<loadfile>>",'"' + comp_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + corel_mod + '"')

	h = open(inlist_corel, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_corel, "inlist")

	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time


def run_reduce(inlist_reduce, corel_mod, reduce_mod, mp):
	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_reduce', 'r')
	g = f.read()
	f.close()


	g = g.replace("<<loadfile>>",'"' + corel_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + reduce_mod + '"')
	g = g.replace("<<mp>>",str((mp * mearth / msun)))
	

	h = open(inlist_reduce, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_reduce, "inlist")


	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time

def run_corem(inlist_corem, reduce_mod, corem_mod, core_mass, rho):
	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_corem', 'r')
	g = f.read()
	f.close()

	g = g.replace("<<loadfile>>",'"' + reduce_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + corem_mod + '"')
	g = g.replace("<<core_mass>>", str(core_mass * mearth / msun))
	g = g.replace("<<rho>>", str(rho))
	
	h = open(inlist_corem, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_corem, "inlist")

	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time


def run_heating(inlist_heating, corem_mod, heating_mod, entropy, luminosity):
	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_heating', 'r')
	g = f.read()
	f.close()

	g = g.replace("<<loadfile>>",'"' + corem_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + heating_mod + '"')
	g = g.replace("<<entropy>>", str(entropy))
	
	# This is to inflate the planet
	g = g.replace("<<luminosity>>", str(luminosity))
	
	h = open(inlist_heating, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_heating, "inlist")

	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time


def run_cooling(inlist_cooling, corem_mod, cooling_mod, entropy):
	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_cooling', 'r')
	g = f.read()
	f.close()

	g = g.replace("<<loadfile>>",'"' + corem_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + cooling_mod + '"')
	g = g.replace("<<entropy>>", str(entropy))
	
	h = open(inlist_cooling, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_cooling, "inlist")

	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time


def run_remove_heating(remove_heating_profile, inlist_remove_heating, heating_mod, remove_mod):
	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_remove_heating', 'r')
	g = f.read()
	f.close()

	g = g.replace("<<loadfile>>",'"' + heating_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + remove_mod + '"')
	g = g.replace("<<remove_heating_profile>>", '"' + remove_heating_profile + '"')
	
	h = open(inlist_remove_heating, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_remove_heating, "inlist")

	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time

def run_remove_cooling(remove_cooling_profile, inlist_remove_cooling, cooling_mod, remove_mod):
	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_remove_cooling', 'r')
	g = f.read()
	f.close()

	g = g.replace("<<loadfile>>",'"' + cooling_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + remove_mod + '"')
	g = g.replace("<<remove_cooling_profile>>", '"' + remove_cooling_profile + '"')
	
	h = open(inlist_remove_cooling, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_remove_cooling, "inlist")

	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time



def run_irrad(irrad_profile, inlist_irrad, remove_mod, irrad_mod, column_depth, flux_dayside):
	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_irrad', 'r')
	g = f.read()
	f.close()

	g = g.replace("<<loadfile>>",'"' + remove_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + irrad_mod + '"')
	g = g.replace("<<irrad_profile>>", '"' + irrad_profile + '"')
	g = g.replace("<<flux_dayside>>", str(flux_dayside))
	g = g.replace("<<column_depth>>", str(column_depth))

	h = open(inlist_irrad, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_irrad, "inlist")

	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time


def run_evolve(evolve_profile, inlist_evolve, irrad_mod, evolve_mod,
				n_frac, a, ms, orb_sep, ec, column_depth, flux_dayside,
				formation_time, teq, BA, escape_regime, diff_sep, homopause_temp):

	start_time = time.time()
	print ("create initial planet")
	f = open('inlist_evolve', 'r')
	g = f.read()
	f.close()

	#File Parameters
	g = g.replace("<<loadfile>>",'"' + irrad_mod + '"')
	g = g.replace("<<smwtfname>>", '"' + evolve_mod + '"')
	g = g.replace("<<evolve_profile>>", '"' + evolve_profile + '"')

	#Flux Parameters
	g = g.replace("<<formation_time>>",str(formation_time))
	g = g.replace("<<column_depth>>",str(column_depth))
	g = g.replace("<<flux_dayside>>", str(flux_dayside))

	#x-controls
	g = g.replace("<<n_frac>>", str(n_frac))
	g = g.replace("<<a>>", str(a))
	g = g.replace("<<ms>>", str(ms))
	g = g.replace("<<BA>>", str(BA))
	g = g.replace("<<orb_sep>>", str(orb_sep))
	g = g.replace("<<ec>>", str(ec))
	g = g.replace("<<teq>>", str(teq))
	g = g.replace("<<escape_regime>>", str(escape_regime))
	g = g.replace("<<diff_sep>>", str(diff_sep))
	g = g.replace("<<homopause_temp>>", str(homopause_temp))
	
	h = open(inlist_evolve, 'w')
	h.write(g)
	h.close()
	shutil.copyfile(inlist_evolve, "inlist")

	os.system('./star_make_planets')
	run_time = time.time() - start_time
	return run_time
