import numpy as np
import matplotlib as mpl
from astropy.table import Table, Column, join 
#mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.wcs import WCS
from astropy.io import fits
import sys
import math
import os
import glob
import sys
from sortedcontainers import SortedDict
import datetime as dt
import imageio
import os
from PIL import Image
from matplotlib.colors import LogNorm
from astropy.nddata.utils import Cutout2D
from astropy import units as u
import datetime as dt 
import astropy.units as u
from astroML.crossmatch import crossmatch_angular 
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, sigma_clip
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import glob
import csv


###-----------------------------------   INPUT FIELD NAME ETC --------------------------------------------------------------###
year = 2016
field =  'ngc6101'
input_path = '/fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/'+str(year)+'/*/'+field+'/g_band/single/*160728*/'
exp_time = 20 
path_list = glob.glob(input_path)

bad_images = []
for i in path_list: 
	try: 
		averages = []
		os.chdir(i)
		#print(i)
		av_path = str(i + 'photom_correction_files/')
		#print(av_path)
		test = os.chdir(av_path)
		#print(os.listdir('.'))
		#print(os.listdir('.'))
		for filename in os.listdir('.'):
			if filename.startswith('average'):
				#print(filename)
				try: 
					avs = np.loadtxt(filename, unpack = True, skiprows= 1)
					#print(avs)
					for f in avs: 
						averages.append(f)
					mean= np.mean(averages)
					#print(mean)
				
					#if mean == np.double: 
				
			
				
				except:
					pass
		a = True
		if a==True:
			#print("YAAS")
			os.chdir("../final_source_cats/")
			for cat in os.listdir('.'):
				if cat.endswith('.cat_NOT_CORRECTED.ascii'):
					
					#print('--------------------------------------------------')
					RA, DEC, g_mag_AUTO, g_mag_err_AUTO, g_mag_APER, g_mag_err_APER = np.loadtxt(cat, unpack=True, skiprows =1)
					zp_exp_correction = 2.5*np.log10(exp_time)
					info_table = Table()
					info_table['RA'] = RA
					info_table['DEC'] = DEC
					info_table['g_mag_AUTO'] = (g_mag_AUTO + zp_exp_correction) - mean
					info_table['g_mag_err_AUTO'] = g_mag_err_AUTO
					info_table['g_mag_APER'] = (g_mag_APER + zp_exp_correction) - mean
					info_table['g_mag_err_APER'] = g_mag_err_APER
			
					t = info_table 
					filename = str(cat)
					cut_filename = filename[:len(filename)-19]
					#print(cut_filename) 
			
					output = cut_filename + 'cat_CORRECTED.ascii'
					print('--------------- NEWLY CORRECTED------------')
					print(output)
					t.write(output, format='ascii', overwrite=True)
						
				
				
				
				
				else: 
					#print('------NOT CORRECTED--------')
					#print(cat)
					print('---------------------------')
				
					pass 
	
		#print(os.listdir('.'))
	
	except:
		pass 
	
bad = Table()
bad['paths'] = bad_images
outputs = '/fred/oz100/NOAO_archive/archive_NOAO_data/scripts/correct_photom/' + str(year) +'bad_images_not_correctable.ascii'	
bad.write(outputs, format='ascii', overwrite=True)		
			
			
	
		#a = os.chdir("../final_source_cats/") 
		#print(a)
		#cats = glob.glob(a)
		#print(str(cats))
			
	#print(averages)
	#print(i)
	#print(np.mean(averages))
