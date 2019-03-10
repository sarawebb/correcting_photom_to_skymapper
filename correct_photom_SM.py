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
import pandas as pd

###-----------------------------------   INPUT FIELD NAME ETC --------------------------------------------------------------

exp_time = 20 
input_path = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/*/SE_cats/'
output_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'
input_cats = '/mnt/dwf/sky_mapperDR2/shortlisted_stars_for_photom/4hr_SM.csv_SHORTLISTED.ascii'

###-----------------------------------^^^^^^^^^^^ DID YOU INPUT FIELD NAMEs ? ^^^^^^ --------------------------------------------------------------


path_list = glob.glob(input_path)
#print(path_list)

for i in path_list: 
	print(i)
	#print(os.listdir(i))
	for filename in os.listdir(i):
		#print(filename)
		if filename.endswith('.cat'):
			## Read in SE catalog from DWF data
			print(i)
			print(filename)
			
			MAG_APER, MAGERR_APER, MAG_AUTO, MAGERR_AUTO, XPEAK_IMAGE, YPEAK_IMAGE, X_IMAGE, Y_IMAGE, ALPHA_J2000, DELTA_J2000 = np.loadtxt(i + filename, unpack = True)
			
			zp_exp_correction = 2.5*np.log10(exp_time)
			
			ccd_SE_cat = Table()
			ccd_SE_cat['mag_aper'] = MAG_APER + zp_exp_correction
			ccd_SE_cat['magerr_aper'] = MAGERR_APER
			ccd_SE_cat['mag_auto'] = MAG_AUTO + zp_exp_correction
			ccd_SE_cat['mag_auto_err'] = MAGERR_AUTO
			ccd_SE_cat['xpeak_image'] = XPEAK_IMAGE
			ccd_SE_cat['ypeak_image'] = YPEAK_IMAGE
			ccd_SE_cat['x_image'] = X_IMAGE
			ccd_SE_cat['y_image'] = Y_IMAGE
			ccd_SE_cat['RA'] = ALPHA_J2000
			ccd_SE_cat['DEC'] = DELTA_J2000	
			
			DWF_X = np.empty((len(ccd_SE_cat), 2), dtype=np.float64)
			DWF_X[:, 0] = ALPHA_J2000
			DWF_X[:, 1] = DELTA_J2000
			#print(ccd_SE_cat)
			
			
			ra, dec, ngood, class_star, g_psf, e_g_psf, r_psf, e_r_psf, i_psf, e_i_psf, z_psf, e_z_psf = np.loadtxt(input_cats , unpack = True, skiprows = 1)
			
			SM_X = np.empty((len(ra), 6), dtype=np.float64)
			SM_X[:, 0] = ra
			SM_X[:, 1] = dec
			SM_X[:, 2] = g_psf
			SM_X[:, 3] = r_psf
			SM_X[:, 4] = i_psf
			SM_X[:, 5] = z_psf
			
			##---------------------------------Cross match DWF CAT with SkyMapper Data ------------------------------- # 
			
			max_radius = 1./3600 #1 arc second 
			
			dist_between, ind_row = crossmatch_angular(DWF_X, SM_X, max_radius)
			match = ~np.isinf(dist_between)
			
			match_table = Table()
			match_table['matched_true_false'] = match 
			match_table['matched_ID'] = ind_row
			match_table['matched_DWF_data_gmag'] = ccd_SE_cat['mag_auto'] = ccd_SE_cat['mag_auto']
			match_table['matched_DWF_RA'] = ccd_SE_cat['RA']
			match_table['matched_DWF_DEC'] = ccd_SE_cat['DEC']
			#print(match_table)
			
			
			SM_match_true = []
			SM_row_matched = []
			DWF_g_mags_matched = []
			DWF_g_mags_error_matched = []
			DWF_obs_ra_matched = []
			DWF_obs_dec_matched = []
				
			for row in match_table: 
				if row['matched_true_false'] == True:
					#print(row)
					SM_match_true.append(row['matched_true_false'])
					SM_row_matched.append(row['matched_ID'])
					DWF_g_mags_matched.append(row['matched_DWF_data_gmag'])
					DWF_obs_ra_matched.append(row['matched_DWF_RA'])
					DWF_obs_dec_matched.append(row['matched_DWF_DEC'])
					
			SM_RA = []
			SM_DEC = []
			SM_g_mag = []
			SM_r_mag = []
			SM_i_mag = []
			SM_z_mag = []
			
			for i in SM_row_matched:
				#print(i)
				
				RA = SM_X[i, 0]
				DEC = SM_X[i, 1]	
				g_mag = SM_X[i, 2]	
				r_mag = SM_X[i, 3]	
				i_mag = SM_X[i, 4]	
				z_mag = SM_X[i, 5]
				SM_RA.append(RA)
				SM_DEC.append(DEC)	
				SM_g_mag.append(g_mag)	
				SM_r_mag.append(r_mag)	
				SM_i_mag.append(i_mag)		
				SM_z_mag.append(z_mag)
					
			
				
			source_match_table = Table()
			source_match_table['SM_match_true'] = SM_match_true
			source_match_table['SM_match_row_ID'] = SM_row_matched
			source_match_table['SM_g_mag'] = SM_g_mag
			source_match_table['SM_r_mag'] = SM_r_mag		
			source_match_table['SM_i_mag'] = SM_i_mag	
			source_match_table['SM_z_mag'] = SM_z_mag	
			source_match_table['SM_RA'] = SM_RA		
			source_match_table['SM_DEC'] = SM_DEC	
			source_match_table['DWF_gmag'] = DWF_g_mags_matched
			source_match_table['DWF_RA'] = DWF_obs_ra_matched		
			source_match_table['DWF_DEC'] = DWF_obs_dec_matched
			#print(source_match_table)
			
			diff_mags = []
			for row, row1 in zip(source_match_table['DWF_gmag'], source_match_table['SM_g_mag']):
				diff = row - row1
				#print(diff)
				check = np.isfinite(diff)
				if check == True:
					diff_mags.append(diff)
				else:
					pass 
			
			#print(diff_mags)
			
			
			mean_diff = np.mean(diff_mags)
			median_diff= np.median(diff_mags)
			
			print('Mean ZP offset: ' + str(mean_diff))
			print('Median ZP offset: ' + str(median_diff))			
			
			
			source_match_table['DWF_g_mag_ZPoffset'] = source_match_table['DWF_gmag'] - mean_diff
			a = [10, 23]
			b = [10, 23]
			plt.scatter(source_match_table['DWF_gmag'], source_match_table['SM_g_mag'])
			plt.ylim(10, 23)
			plt.xlim(10, 23)
			plt.plot(a, b, 'r')
			plt.xlabel('DWF_g_mag')
			plt.ylabel('SM_g_mag')
			plt.title('uncorrected')
			plt.show()
			
			a = [10, 23]	
			b = [10, 23]
			plt.scatter(source_match_table['DWF_g_mag_ZPoffset'], source_match_table['SM_g_mag'])
			plt.ylim(10, 23)
			plt.xlim(10, 23)
			plt.plot(a, b, 'r')
			plt.xlabel('DWF_g_mag')
			plt.ylabel('SM_g_mag')
			plt.title('corrected')
			plt.show()
			
			
			
