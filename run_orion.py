#!/usr/bin/env python
#python script created by Phil Arras UVA, modified by HW Chen NU and then modified more by Isaac Malsky
import numpy as np
import os
import scipy
from scipy import loadtxt, optimize
import mysubsprograms as my
import sys

#Constants
msun = 1.9892e33
rsun = 6.9598e10
rearth = 6.371008e8
mjup = 1.8986e30
rjup = 6.9911e9
mearth = 5.97e27
sigma=5.67e-5
au = 1.496e13

initial_mod = "initial_planet.mod"
####################################################
#########         PARAMETERS LISTS         #########
####################################################
mpList=[6.261]
orbitalList=[0.01411]
enFracList=[0.012]
yList = [.24]
zList = [.02]
entropyList = [-1]
n_frac_list = [0.10]

####################################################
#########        Star Parameters              ######
####################################################                            
rs = 0.216                      #star radius in rsun
ms = 0.150                           #host_star_mass
Teff_star = 3026                #Host Star Temp

####################################################
#########        hompause temperature         ######
####################################################   
homopause_temp = 10000


####################################################
#########     Mass Loss Assumptions           ######
####################################################    
a = 1.0                          #frac_absorbing_radius
BA= 0.20                         #planet Bond albedo
ec = 1e9                         #eddy coefficient
formation_time = 6e6             #Disk formation time

# Mass Loss Regime
escape_type = 'hu'

# Diffusive Separation, 1 is on 0 is off
diff_sep = 1

# This is for setting the type of escape
# Only hu et. at 2017 or something
# And Murray Clay et. al or something are programmed in
if escape_type == 'hu':
    escape_regime = 0
elif escape_type == 'murray':
    escape_regime = 1
else:
    escape_regime = None


####################################################
#########                 Run!                ######
####################################################  
for mp in mpList:
    pre_reduce_mod = "pre_reduce_" + str(mp) + "_" + ".mod"
    inlist_pre_reduce = "inlist_pre_reduce_" + str(mp)
    run_time = my.run_pre_reduce(inlist_pre_reduce, initial_mod, pre_reduce_mod, mp)

    for enFrac in enFracList:
        core_mass = mp * (1 - enFrac)
        rho = my.calculate_rho(mp, enFrac)[0]
        radius1 = my.calculate_rho(mp, enFrac)[1]

        pre_core_mod = "pre_core_" + str(mp) + "_" + str(enFrac) + ".mod"
        inlist_pre_core = "inlist_pre_core_" + str(mp) + "_" + str(enFrac)
        run_time = my.run_pre_core(inlist_pre_core, pre_reduce_mod, pre_core_mod, enFrac,core_mass,rho)

        for z in zList:
            for y in yList:
                comp_mod = "comp_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + ".mod"
                inlist_comp = "inlist_comp_" + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z)
                run_time = my.run_comp(inlist_comp, pre_core_mod, comp_mod, z, y)

                corel_mod = "corel_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + ".mod"
                inlist_corel = "inlist_corel_"  + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z)
                run_time = my.run_corel(inlist_corel, comp_mod, corel_mod)

                reduce_mod = "reduce_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + ".mod"
                inlist_reduce = "inlist_reduce_" + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z)
                run_time = my.run_reduce(inlist_reduce, corel_mod, reduce_mod, mp)

                corem_mod = "corem_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + ".mod"
                inlist_corem = "inlist_corem_" + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z)
                run_time = my.run_corem(inlist_corem, reduce_mod, corem_mod, core_mass, rho)


                for entropy in entropyList:
                    if (os.path.isfile('LOGS/' + corem_mod) == True):
                        entropy_list, luminosity_list = loadtxt('LOGS/' + corem_mod, unpack=True, skiprows =6, usecols=[12,13])
                        
                        if isinstance(entropy_list, (list, np.ndarray)):
                            currentropy = entropy_list[-1]
                            luminosity = luminosity_list[-1]
                        else:
                            currentropy = entropy_list
                            luminosity = luminosity_list
                    else:
                        break

                    # Try a couple ways of fitting entropy
                    if entropy == -1:
                        entropy = np.round(7.0 + (mp / 25.0), 2)
                    else:
                        entropy = entropy

                    # THIS NEEDS TO BE MODIFIED TO CORRECTLY INFLATE PLANTES
                    #luminosity = ((0.01 * (mp ** 2.0)) + (0.10 * mp) + 25) * luminosity
                    luminosity = 25 * luminosity

                    if currentropy < float(entropy):
                        heating_mod = "heating_" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy) + ".mod"
                        inlist_heating = "inlist_heating_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(entropy)
                        run_time = my.run_heating(inlist_heating, corem_mod, heating_mod, entropy, luminosity)

                        remove_mod = "remove_heating_" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy) + ".mod"
                        remove_heating_profile = "profile_remove_heating" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy)
                        inlist_remove_heating = "inlist_remove_heating_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(entropy)
                        run_time = my.run_remove_heating(remove_heating_profile, inlist_remove_heating, heating_mod, remove_mod)


                        for orb_sep in orbitalList:
                            for n_frac in n_frac_list:
                                if (os.path.isfile(remove_mod) == True):

                                    flux_dayside = (sigma*Teff_star**4 * (rs * rsun / orb_sep / au )**2)*(1-BA)    # flux hitting planet's dayside
                                    teq = (flux_dayside / 4.0 / sigma) ** 0.25    #equalibrium temperature
                                    column_depth = my.calculate_column_depth(teq, remove_heating_profile, Teff_star)

                                    if (column_depth != -1):
                                        irrad_mod = "irrad_" + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy)+ ".mod"
                                        inlist_irrad = "inlist_irrad_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)
                                        irrad_profile = "profile_irrad" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac)
                                        run_time = my.run_irrad(irrad_profile, inlist_irrad, remove_mod, irrad_mod, column_depth, flux_dayside)
                                        
                                        evolve_mod = "evolve_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac) + ".mod"
                                        evolve_profile = "profile_evolve" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac)
                                        inlist_evolve = "inlist_evolve_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep) + "_" + str(entropy) + "_" + str(n_frac)
                                        run_time = my.run_evolve(evolve_profile, inlist_evolve, irrad_mod,
                                            evolve_mod, n_frac, a, ms, orb_sep, ec, column_depth, flux_dayside,
                                            formation_time, teq, BA, escape_regime, diff_sep, homopause_temp)
                                    else:
                                        pass
                                else:
                                    pass


                    else:
                        cooling_mod = "cooling_" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy) + ".mod"
                        inlist_cooling = "inlist_cooling_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(entropy)
                        run_time = my.run_cooling(inlist_cooling, corem_mod, cooling_mod, entropy)

                        remove_mod = "remove_cooling_" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy) + ".mod"
                        remove_cooling_profile = "profile_remove_cooling_" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy)
                        inlist_remove_cooling = "inlist_remove_cooling_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(entropy)
                        run_time = my.run_remove_cooling(remove_cooling_profile, inlist_remove_cooling, cooling_mod, remove_mod)


                        for orb_sep in orbitalList:
                            for n_frac in n_frac_list:
                                if (os.path.isfile(remove_mod) == True):
                                    flux_dayside = (sigma * Teff_star ** 4 * (rs * rsun / orb_sep / au ) ** 2) * (1 - BA)    # flux hitting planet's dayside
                                    teq = (flux_dayside/ 4.0 /sigma)**0.25    #equalibrium temperature
                                    column_depth = my.calculate_column_depth(teq, remove_cooling_profile, Teff_star)

                                    if (column_depth != -1):
                                        irrad_mod = "irrad_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + ".mod"
                                        inlist_irrad = "inlist_irrad_" + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z) + "_" + str(orb_sep) + "_" + str(entropy)
                                        irrad_profile = "profile_irrad" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac)
                                        run_time = my.run_irrad(irrad_profile, inlist_irrad, remove_mod, irrad_mod, column_depth, flux_dayside)

                                        evolve_mod = "evolve_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac) + ".mod"
                                        evolve_profile = "profile_evolve" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac)
                                        inlist_evolve = "inlist_evolve_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep) + "_" + str(entropy) + "_" + str(n_frac)

                                        run_time = my.run_evolve(evolve_profile, inlist_evolve, irrad_mod,
                                            evolve_mod, n_frac, a, ms, orb_sep, ec, column_depth, flux_dayside,
                                            formation_time, teq, BA, escape_regime, diff_sep, homopause_temp)
                                    else:
                                        pass
                                else:
                                    pass
