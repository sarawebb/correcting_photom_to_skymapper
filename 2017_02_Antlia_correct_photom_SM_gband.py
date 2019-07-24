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


###-----------------------------------   INPUT FIELD NAME ETC --------------------------------------------------------------

exp_time = 20 
input_path = '/mnt/dwf/archive_NOAO_data/data_outputs/2017/02/Antlia/g_band/single/*/SE_cats/'
input_cats = '/mnt/dwf/sky_mapperDR2/shortlisted_stars_for_photom/Antlia_SkyM.csv_SHORTLISTED.ascii'

###-----------------------------------^^^^^^^^^^^ DID YOU INPUT FIELD NAMEs ?  --------------------------------------------------------------




path_list = glob.glob(input_path)

for i in path_list: 
	print(i)
	av_correction = []
	odd_differences = []
	filenames_good = [] 
	filenames_odd = [] 
	for filename in os.listdir(i):
		#av_correction = np.zeros(len(os.listdir(i)))
		#print(av_correction)
		#print( filename)
		#filenames_for_av = []
		#filenames_for_av.append(filename)
		if filename.endswith('.cat'):
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
			ra, dec, ngood, class_star, g_psf, e_g_psf, r_psf, e_r_psf, i_psf, e_i_psf, z_psf, e_z_psf = np.loadtxt(input_cats , unpack = True, skiprows = 1)
			SM_X = np.empty((len(ra), 6), dtype=np.float64)
			SM_X[:, 0] = ra
			SM_X[:, 1] = dec
			SM_X[:, 2] = g_psf
			SM_X[:, 3] = r_psf
			SM_X[:, 4] = i_psf
			SM_X[:, 5] = z_psf
			os.chdir(i)
			os.chdir("..")
			filename_directory = os.path.abspath(os.curdir)
			photom_correction_path = filename_directory + '/photom_correction_files/'
			#print(photom_correction_path)
			final_source_cat_path = filename_directory + '/final_source_cats/'
			#print(final_source_cat_path)
			photom_correction_path_images = photom_correction_path + 'skymapper_checkphotom_plots/'
			#print(photom_correction_path_images)
			
			if not os.path.exists(photom_correction_path):
		               	os.makedirs(photom_correction_path, 0o755)
			else:
				pass 
			if not os.path.exists(final_source_cat_path):
		               	os.makedirs(final_source_cat_path, 0o755)
			else:
		               	pass
			if not os.path.exists(photom_correction_path_images):
		               	os.makedirs(photom_correction_path_images, 0o755)
			else:
				pass 
			
			max_radius = 1./3600 #1 arc second 
			dist_between, ind_row = crossmatch_angular(DWF_X, SM_X, max_radius)
			match = ~np.isinf(dist_between)
			if len(match) != 0: 
				match_table = Table()
				match_table['matched_true_false'] = match
				match_table['matched_ID'] = ind_row
				match_table['matched_DWF_data_gmag'] = ccd_SE_cat['mag_auto']
				match_table['matched_DWF_RA'] = ccd_SE_cat['RA']
				match_table['matched_DWF_DEC'] = ccd_SE_cat['DEC']
				SM_match_true = []
				SM_row_matched = []
				DWF_g_mags_matched = []
				DWF_g_mags_error_matched = []
				DWF_obs_ra_matched = []
				DWF_obs_dec_matched = []
				for row in match_table:
					if row['matched_true_false'] == True:
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
				for j in SM_row_matched:
					RA = SM_X[j, 0]
					DEC = SM_X[j, 1]
					g_mag = SM_X[j, 2]
					r_mag = SM_X[j, 3]
					i_mag = SM_X[j, 4]
					z_mag = SM_X[j, 5]
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
				diff_mags = []
				for row, row1 in zip(source_match_table['DWF_gmag'], source_match_table['SM_g_mag']):
					if row1 < 24:
						diff = row - row1
						if diff < 1.5: 
							diff_mags.append(diff)
							filenames_good.append(filename)
							#print(diff)
						elif diff > 1.5: 
							odd_differences.append(diff)
							filenames_odd.append(filename)
					elif row1 > 24:
						pass
				if len(diff_mags) > 0: 
					mean_diff = np.mean(diff_mags)
					median_diff= np.median(diff_mags)
					#print('got to step 3: found mean: ' + str(mean_diff))
					source_match_table['DWF_g_mag_ZPoffset'] = source_match_table['DWF_gmag'] - mean_diff
					a = [10, 23]
					b = [10, 23]
					plt.scatter(source_match_table['DWF_gmag'], source_match_table['SM_g_mag'])
					plt.ylim(10, 23)
					plt.xlim(10, 23)
					plt.plot(a, b, 'r')
					plt.xlabel('DWF_g_mag')
					plt.ylabel('SM_g_mag')
					plt.title('uncorrected ' + str(filename))
					plt.savefig(photom_correction_path_images + filename + 'skymapper_uncorrection.png', overwrite = True)
					#plt.show()
					plt.close()
					plt.scatter(source_match_table['DWF_g_mag_ZPoffset'], source_match_table['SM_g_mag'])
					plt.ylim(10, 23)
					plt.xlim(10, 23)
					plt.plot(a, b, 'r')
					plt.xlabel('DWF_g_mag')
					plt.ylabel('SM_g_mag')
					plt.title('corrected ' + str(filename))
					plt.savefig(photom_correction_path_images + filename + 'skymapper_corrected.png', overwrite = True)
					plt.close()
					
					corrected_g_mags = ccd_SE_cat['mag_auto'] - mean_diff
					corrected_g_mags_aper = ccd_SE_cat['mag_aper']  - mean_diff
					#print(corrected_g_mags)
					#print('Matching Sky Mapper sources found! Correcting to average deviation for this ccd.')
					new_cat = Table() 
					new_cat['RA'] = ccd_SE_cat['RA']
					new_cat['DEC'] = ccd_SE_cat['DEC']
					new_cat['g_mag'] = corrected_g_mags
					new_cat['g_mag_err'] = ccd_SE_cat['mag_auto_err']
					new_cat['g_mag_APER'] = corrected_g_mags_aper
					new_cat['g_mag_err_APER'] = ccd_SE_cat['magerr_aper']
					#print(new_cat)
					av_correction.append(mean_diff)
					#print('HELOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO')
					#print('average correction ' + str(av_correction))
					t = new_cat
					output= final_source_cat_path + '/'+ filename + '_CORRECTED.ascii'
					#print(output)
					t.write(output, format= 'ascii' , overwrite = True)
				else: 
					print('No gband from SkyMapper for this ccd')
					if len(av_correction) != 0:
						mean_of_means = np.mean(av_correction)
						print('mean_of_means: ' + str(mean_of_means))
					
						if mean_of_means == 0: 
							print('No matching Skymapper sources, NOT correcting for this ccd')
							corrected_g_mags = ccd_SE_cat['mag_auto']  - mean_of_means
							corrected_g_mags_aper = ccd_SE_cat['mag_aper']  - mean_of_means
							new_cat = Table()
							new_cat['RA'] = ccd_SE_cat['RA']
							new_cat['DEC'] = ccd_SE_cat['DEC']
							new_cat['g_mag_AUTO'] = corrected_g_mags
							new_cat['g_mag_err_AUTO'] = ccd_SE_cat['mag_auto_err']
							new_cat['g_mag_APER'] = corrected_g_mags_aper
							new_cat['g_mag_err_APER'] = ccd_SE_cat['magerr_aper']
						
							t = new_cat
							output= final_source_cat_path + '/' + filename + '_NOT_CORRECTED.ascii'
							t.write(output, format= 'ascii' , overwrite = True)
						
						elif mean_of_means > 0:
							#print('No matching Skymapper sources, correcting to average of average corrections for this field')
							corrected_g_mags = ccd_SE_cat['mag_auto']  - mean_of_means
							corrected_g_mags_aper = ccd_SE_cat['mag_aper']  - mean_of_means
							new_cat = Table() 
							new_cat['RA'] = ccd_SE_cat['RA']
							new_cat['DEC'] = ccd_SE_cat['DEC']
							new_cat['g_mag_AUTO'] = corrected_g_mags
							new_cat['g_mag_err_AUTO'] = ccd_SE_cat['mag_auto_err']
							new_cat['g_mag_APER'] = corrected_g_mags_aper
							new_cat['g_mag_err_APER'] = ccd_SE_cat['magerr_aper']
							t = new_cat
							output= final_source_cat_path + '/' + filename + '_CORRECTED.ascii'
							t.write(output, format= 'ascii' , overwrite = True)
					else: 
						print('No matching Skymapper sources, NOT correcting for this ccd')
						corrected_g_mags = ccd_SE_cat['mag_auto']  
						corrected_g_mags_aper = ccd_SE_cat['mag_aper']
						new_cat = Table()
						new_cat['RA'] = ccd_SE_cat['RA']
						new_cat['DEC'] = ccd_SE_cat['DEC']
						new_cat['g_mag_AUTO'] = corrected_g_mags
						new_cat['g_mag_err_AUTO'] = ccd_SE_cat['mag_auto_err']
						new_cat['g_mag_APER'] = corrected_g_mags_aper
						new_cat['g_mag_err_APER'] = ccd_SE_cat['magerr_aper']

						t = new_cat
						output= final_source_cat_path + '/' + filename + '_NOT_CORRECTED.ascii'
						t.write(output, format= 'ascii' , overwrite = True)
					
		               	
		               	
		        
		
	try: 		
		av_tab = Table()
		av_tab['average correction'] = av_correction
		
		#av_tab['filenames'] = filenames_good
		print('HELLOOOOOOOOOOOOOOOOOOOOOOOOOOOO')
		print(av_tab)
		#av_tab['filename'] = filenames_for_av
			
		output2 = photom_correction_path + '/average_correction_per_ccd.ascii'
		av_tab.write(output2, format ='ascii', overwrite = True) 
		#print(av_tab)
	except: 
		pass 
		
	try: 	
		odd_tab = Table()
		odd_tab['odd differences'] = odd_differences
		odd_tab['filenames'] = filenames_odd
		output3 = photom_correction_path + '/odd_values_found_per_ccd.ascii'
		odd_tab.write(output3, format ='ascii', overwrite = True)
	except: 
		pass 		
		
				
	filenames_good.clear()		
	av_correction.clear()			
			
			
			
				
			
