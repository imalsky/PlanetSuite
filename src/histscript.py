#!/usr/bin/env python

import math
import numpy as np
import os
import shutil
import mysubsprograms as my
import sys
import random
from math import log10

msun = 1.9892e33
rsun = 6.9598e10
mjup = 1.8986e30
rjup = 6.9911e9
mearth = 5.97e27
sigma=5.67e-5
au = 1.496e13

# fiducial parameters
mp = 0.055  				# total planet mass in mjup
rp = 3.0				# inital planet radius in rjup
mcore = 10.0				# planet core mass in mearth
rhocore = 10.0				# core density in g/cc
mp_wo_core = mp - mcore*mearth/mjup	# mass of planet before core put in
z = 0.03                                # metallicity of both planet and star
y = 0.25                                # helium fraction of both planet and (initial) star
maxage = 1.e10                          # ending age
initialEntropy= 8                    # entropy in kb/Nb, after the planet formation
coolingTime = 1e6          # cooling timescale before mass loss evolution

# parameters having to do with irradiation
ms = 1.0				# star mass in msun
rs = 1.0				# star radius in rsun
Teff_star = 5800			# stellar Teff
orb_sep = 0.1   # orbital sepration in AU (to set day-side flux)  134000,13400000,1340000000 ergs/s/cm^2   ie. 0.032, 0.32,3.2
BA= 0.1
irrad_col = 300.0			# column depth for depositing stellar radiation as heat
flux_dayside = (sigma*Teff_star**4 * ( rs*rsun / orb_sep / au )**2)*(1-BA)	# flux hitting planet's dayside
Teq = (flux_dayside/4.0/sigma)**0.25
Mmin=1                                  #for the random generator
Mmax=100
enFracmin=0.001
enFracmax=0.99
minitial=30
# flags to skip steps
do_create_planet = True
do_put_in_core = True
do_set_entropy= True
do_cool = True
do_relaxm = True
do_evolve_planet = True
do_evolve_planet2 = True

f=open('logfile','w')

mpList=[5.7,7,10,15,20,25,30,35,40,45]
enFracList=[0.22,0.001,0.002,0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.25]

mpList=[4,5,6,7,8,9,11,12,13,15,17,20,24,26,30]

enFracList=[0.25,0.4,0.5,0.571,0.625,0.6666,0.7273,0.75,0.769,0.8,0.824,0.85,0.875,0.885,0.9]



mC=[]
rC=[]
TempEq=[]
Kappa=[]
oneBar=1E+05
logMList=[]
logenFracList=[]
MList=[]
enFracList=[]

for i in range(1200):
	num=random.random()
	logM=log10(Mmin) + num*(log10(Mmax) - log10(Mmin))
	logMList.append(logM)
 
mpList=[10**x for x in logMList]

for i in range(1200):
	num=random.random()
	logenFrac=log10(enFracmin) + num*(log10(enFracmax) - log10(enFracmin))
	logenFracList.append(logenFrac)

enFracList=[10**x for x in logenFracList]





#for i in range(0,40):
mpList=[10,100]

for i in range(0,1200):
	if mpList[i] <= 30:
		createmodel = "planet_create_" + str(minitial) + "_ME_" + str(rp) + "_RJ.mod"

		inlist1 = "inlist_create_" + str(mp)[0:6] + "_MJ_" + str(rp) + "_RJ"
		run_time = my.create_planet(mp,y,z,inlist1,createmodel)


		k=open('LOGS/history.data','r')
		for line in k.readlines():
			pass
		last_temp=line
		last=last_temp.split()

                Rmp=mpList[i]
		if do_put_in_core:
			enFrac=enFracList[i]
			mcore=mpList[i]-mpList[i]*enFrac
             
			my.print_parameters(Rmp,enFrac,rhocore,mcore,irrad_col,flux_dayside,Teq,y,z,maxage)
			with open("coreMRcomp2_v40_all.txt",'r') as f:   #determining core radius base on mass, from LA Rogers's model   
				for l in range(11):
					next(f) 
				for line in f:
					mC.append(line.split(None, 10)[0])
					rC.append(line.split(None, 10)[1])
				for k in range(len(mC)):
					mC[k] = float(mC[k]) 
				for k in range(len(rC)):
					rC[k] = float(rC[k])
		
			from scipy.interpolate import interp1d

			#f2=interp1d(mC,rC)
			try:
				f2=interp1d(mC,rC)
				coreRadius= f2(mcore)*6.37e8
				rhocore= (mcore*5.97e27)/((4/3)*(6.37e8*f2(mcore))**3*3.14)
			except:
				do_relaxm = False
				do_evolve_planet = False
				do_evolve_planet2 = False
				
			

			coremodel = "planet_core_" + str(Rmp)[0:6] + "_ME_" + str(mcore)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"
			relaxedmod = "planet_relaxedM_" + str(Rmp)[0:6] + "_ME_" + str(mcore)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"
                        entropymodel = "planet_setEntropy_" + str(Rmp)[0:6] + "_ME_" + str(mcore)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"
			entropymod = "planet_setEntropy_" + str(Rmp)[0:6] + "_ME_" + str(mcore)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"
			coolingmod = "planet_coolingTime_" + str(Rmp)[0:6] + "_ME_" + str(mcore)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"
			evolvemod = "planet_evolve_" + str(Rmp)[0:6] + "_ME_" + str(mcore)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"
                        irradmodel = "planet_removeC_" + str(Rmp)[0:6] + "_ME_" + str(mcore)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"
			evolvemodel2 = "planet_evolve2_" + str(Rmp)[0:6] + "_ME_" + str(mcore)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"

			inlist2 = "inlist_core_" + str(Rmp)[0:6] + "_ME_" + str(enFrac)[0:6]  
			run_time = my.put_core_in_planet(mcore,rhocore,coreRadius,inlist2,createmodel,coremodel)

		if not os.path.exists(coremodel):
			continue


               

        	if do_relaxm:
			inlist3 = "inlist_relaxm_" + str(Rmp)[0:6] + "_ME_" + str(enFrac)[0:6] + "_ME_"
			run_time = my.relaxm(Rmp,inlist3,coremodel,relaxedmod)

		with open('FreedmanEt2014', 'r') as inF:
          		for line in inF:
				line1=line.split(None,14)[1]
            			if float(oneBar) == float(line1) :
					TempEq.append(line.split(None, 14)[0])
					Kappa.append(line.split(None, 14)[3])
			for k in range(len(TempEq)):
				TempEq[k] = float(TempEq[k]) 
			for k in range(len(Kappa)):
			        Kappa[k] = float(Kappa[k])
			#f3=interp1d(TempEq,Kappa)
			#irrad_col= 1/f3(Teq)

                knob= ".false."

                maxage= 6e7


		if do_evolve_planet:
			inlist5 = "inlist_evolve_"  + str(Rmp)[0:6] + "_ME_" + str(enFrac)[0:6]  
			run_time = my.evolve_planet(Teq,irrad_col,flux_dayside,maxage,inlist5,relaxedmod,evolvemod,orb_sep,Rmp,enFrac,knob)

		knob= ".true."
		maxage= 6e9

		if do_evolve_planet2:
			inlist5 = "inlist_evolve_"  + str(Rmp)[0:6] + "_ME_" + str(enFrac)[0:6]  
			run_time = my.evolve_planet(Teq,irrad_col,flux_dayside,maxage,inlist5,evolvemod,evolvemodel2,orb_sep,Rmp,enFrac,knob)




#######below for super-Neptune mass simulations######
	else:     
     		mp=mpList[i]
		enFrac=enFracList[i]
		mcore= 50  #mpList[i]-mpList[i]*enFrac

		my.print_parameters(mp,enFrac,rhocore,mcore,irrad_col,flux_dayside,Teq,y,z,maxage)

		createmodel = "planet_create_" + str(mp)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"

		inlist1 = "inlist_create_" + str(mp)[0:6] + "_MJ_" + str(rp) + "_RJ"
		run_time = my.create_planet(mp,y,z,inlist1,createmodel)

		success = True
		if not os.path.exists(createmodel):
			success=False	
		k=open('LOGS/history.data','r')
		for line in k.readlines():
		  pass
		last_temp=line
		last=last_temp.split()
		print "final model number in create=",last[0]
		print "last[0]==1000",last[0]=="1000"
		if last[0]=="1000":
	 	 success=False
		outstring = '%6.3f\t%6.3f\t%6.3f\t%s\n' % ( mp, rp, run_time, success )

        


		with open("coreMRcomp3_v40_all.txt",'r') as f:   #determining core radius base on mass, from LA Rogers's model     
			for l in range(11):
				next(f) 
			for line in f:
				mC.append(line.split(None, 11)[0])
				rC.append(line.split(None, 11)[6])
			for k in range(len(mC)):
				mC[k] = float(mC[k]) 
			for k in range(len(rC)):
				rC[k] = float(rC[k])

		from scipy.interpolate import interp1d
		f2=interp1d(mC,rC)
		coreRadius= f2(mcore)*6.37e8
		rhocore= (mcore*5.97e27)/((4/3)*(6.37e8*f2(mcore))**3*3.14)

		coremodel = "planet_core_" + str(mp)[0:6] + "_MJ_" + str(mcore)[0:6] + "_ME_" + str(rp) + "_RJ.mod"
		relaxedmod = "planet_relaxedM_" + str(mp)[0:6] + "_ME_" + str(mcore)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"
		entropymodel = "planet_setEntropy_" + str(mp)[0:6] + "_MJ_" + str(mcore)[0:6] + "_ME_" + str(rp) + "_RJ.mod"
 		coolingmodel = "planet_coolingTime_" + str(mp)[0:6] + "_MJ_" + str(mcore)[0:6] + "_ME_" + str(rp) + "_RJ.mod"
		irradmodel = "planet_irrad_" + str(mp)[0:6] + "_MJ_" + str(mcore)[0:6] + "_ME_" + str(rp) + "_RJ.mod"
		evolvemodel2 = "planet_evolve2_" + str(mp)[0:6] + "_MJ_" + str(mcore)[0:6] + "_ME_" + str(rp) + "_RJ.mod"
		evolvemod = "planet_evolve_" + str(mp)[0:6] + "_ME_" + str(mcore)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"

		inlist2 = "inlist_core_" + str(mp)[0:6] + "_MJ_" + str(mcore)[0:6] + "_ME_" + str(rp) + "_RJ"
		run_time = my.put_core_in_planet(mcore,rhocore,coreRadius,inlist2,createmodel,coremodel)



        	if do_relaxm:
			inlist3 = "inlist_relaxm_" + str(mp)[0:6] + "_ME_" + str(enFrac)[0:6] + "_ME_"
			run_time = my.relaxm(mp,inlist3,coremodel,relaxedmod)

		with open('FreedmanEt2014', 'r') as inF:
          		for line in inF:
				line1=line.split(None,14)[1]
            			if float(oneBar) == float(line1) :
					TempEq.append(line.split(None, 14)[0])
					Kappa.append(line.split(None, 14)[3])
			for k in range(len(TempEq)):
				TempEq[k] = float(TempEq[k]) 
			for k in range(len(Kappa)):
			        Kappa[k] = float(Kappa[k])
			#f3=interp1d(TempEq,Kappa)
			#irrad_col= 1/f3(Teq)

                knob= ".false."

                maxage= 6e7


		if do_evolve_planet:
			inlist5 = "inlist_evolve_"  + str(mp)[0:6] + "_ME_" + str(enFrac)[0:6]  
			run_time = my.evolve_planet(Teq,irrad_col,flux_dayside,maxage,inlist5,coremodel,evolvemod,orb_sep,mp,enFrac,knob)

		knob= ".true."
		maxage= 6e9

		if do_evolve_planet2:
			inlist5 = "inlist_evolve_"  + str(mp)[0:6] + "_ME_" + str(enFrac)[0:6]  
			run_time = my.evolve_planet(Teq,irrad_col,flux_dayside,maxage,inlist5,evolvemod,evolvemodel2,orb_sep,mp,enFrac,knob)


f.close()
