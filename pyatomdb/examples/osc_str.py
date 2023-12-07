## Author: Priyanka Chakraborty 07/30/2022
## Description: this code add an extra column for oscillator strength to apec_v3.0.9_line.fits


import numpy as np
from astropy.io import fits
import pandas as pd
import time
from itertools import compress
import math
#import pyatomdb
import os
#import atomdb,apec

t = time.process_time()



if os.path.exists(os.path.expandvars('$ATOMDB/apec_v3.0.9_osc.fits')) == False:

	emi_file = fits.open(os.path.expandvars('$ATOMDB/apec_v3.0.9_line.fits'))  ## opening fits file

	for j in range(2,len(emi_file)):  ### running loops on emissivity extension ##

		emi_lambda = emi_file[j].data['Lambda']
		emi_elem = emi_file[j].data['Element']
		emi_ion = emi_file[j].data['Ion']
		emi_ulev = emi_file[j].data['UpperLev']
		emi_llev = emi_file[j].data['LowerLev']

		
		fij = []
		datacache={}

		settings=apec.parse_par_file(os.path.expandvars('$ATOMDB/apec.par'))


		for i in range (len(emi_elem)):
			index = 0.1
			Z = emi_elem[i]
			z1 = emi_ion[i]

			### fetching upperlevel and lower level from ladat and lvdat
			ladat = atomdb.get_data(Z,z1,'LA', datacache=datacache, settings=settings)
			lvdat = atomdb.get_data(Z,z1,'LV', datacache=datacache, settings=settings)

			if ladat != False:
			
				upper_lev = ladat[1].data['upper_lev']
				lower_lev = ladat[1].data['lower_lev']
			
				#### fetiching indices for matched upper and lower level
				get_indecies = lambda xx, xxs: [iii for (yy, iii) in zip(xxs, range(len(xxs))) if xx == yy]
			
				index_up = get_indecies(emi_ulev[i],upper_lev)
				index_lo = get_indecies(emi_llev[i],lower_lev)
			
				for ii in range(len(index_up)):
					if index_up[ii] in index_lo:
						index = index_up[ii]
						break
					
				if index != 0.1:
				
					#### estimating oscillator strength ####
					larate = ladat[1].data['EINSTEIN_A'][index]
					wavelength = ladat[1].data['wavelen'][index]
		
					stat_weight = lvdat[1].data['lev_deg'][upper_lev[index]-1]/(lvdat[1].data['lev_deg'][lower_lev[index]-1])
					fij.append(1.4992e-8*larate*1e-8*(wavelength**2)*stat_weight)
			
				else:
					larate = 0.0
					wavelength = 0.0
		
					stat_weight = 0.0
					fij.append(0.0)
				
				print(j,i,Z,z1,emi_ulev[i],emi_llev[i])

		#### writing oscillator strength to fits file as a new column ####    
		osc_st = fits.Column(name='Oscil_str', array=fij,format='1E')
		new_column = emi_file[j].columns + osc_st
		new_hdu = fits.BinTableHDU.from_columns(new_column)
		emi_file[j].data = new_hdu.data
		emi_file.writeto('apec_v3.0.9_osc.fits',overwrite=True)

		print('elapsed time for',j,'th extension',time.process_time() - t)



	print('time for entire run',time.process_time() - t)

	emi_file.close()

else: 

	## limiting upper wavelength limit to 40 angstorm

	linefile = "$ATOMDB/apec_v3.0.9_osc.fits"

	lfile = os.path.expandvars(linefile)
	linedata = fits.open(lfile)

	linedata[1].header['COMMENTS'] = 'The upper-limit of wavelength is 40 angstrom'

	ext_num = len(linedata[1].data['kT']) ## number of temperature parameters

	for i in range(2, ext_num):

		wlambda = linedata[i].data['Lambda']
		wlambda_err = linedata[i].data['Lambda_Err']
		epsilon = linedata[i].data['Epsilon']
		epsilon_err = linedata[i].data['Epsilon_Err']
		element = linedata[i].data['Element']
		ion = linedata[i].data['Ion']
		upperlev = linedata[i].data['UpperLev']
		lowerlev = linedata[i].data['LowerLev']
		oscil_str = linedata[i].data['Oscil_str']

		## getting the index where lambda > 40 ang

		index = np.where(wlambda < 40)[0]

		wlambda_new = fits.Column(name='Lambda', array=wlambda[index],format='1E')
		wlambda_err_new = fits.Column(name='Lambda_Err', array=wlambda_err[index],format='1E')
		epsilon_new = fits.Column(name='Epsilon', array=epsilon[index],format='1E')
		epsilon_err_new = fits.Column(name='Epsilon_Err', array=epsilon_err[index],format='1E')
		element_new = fits.Column(name='Element', array=element[index],format='1J')
		ion_new = fits.Column(name='Ion', array=ion[index],format='1J')
		upperlev_new = fits.Column(name='UpperLev', array=upperlev[index],format='1J')
		lowerlev_new = fits.Column(name='LowerLev', array=lowerlev[index],format='1J')
		oscil_str_new = fits.Column(name='Oscil_str', array=oscil_str[index],format='1E')

		new_hdu = fits.BinTableHDU.from_columns([wlambda_new, wlambda_err_new, epsilon_new, epsilon_err_new, element_new, ion_new, upperlev_new, lowerlev_new, oscil_str_new])
		linedata[i].data = new_hdu.data
		


	linedata.writeto(os.path.expandvars('$ATOMDB/apec_v3.0.9_osc.fits'), overwrite=True)






