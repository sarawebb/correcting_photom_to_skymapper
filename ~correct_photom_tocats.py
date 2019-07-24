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

input_path = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/c4d_150116_053736_ooi_g_v1/SE_cats/'
output_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'
output_cats = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/4hr/vizier_cats/'

path_list = glob.glob(input_path)


for i in path_list: 
	for filename in os.listdir(i):
		if filename.endswith('.cat'):
			## Read in SE catalog from DWF data
			MAG_APER, MAGERR_APER, MAG_AUTO, MAGERR_AUTO, XPEAK_IMAGE, YPEAK_IMAGE, X_IMAGE, Y_IMAGE, ALPHA_J2000, DELTA_J2000 = np.loadtxt(i + filename, unpack = True)
			
			ccd_SE_cat = Table()
			ccd_SE_cat['mag_aper'] = MAG_APER
			ccd_SE_cat['magerr_aper'] = MAGERR_APER
			ccd_SE_cat['mag_auto'] = MAG_AUTO
			ccd_SE_cat['mag_auto_err'] = MAGERR_AUTO
			ccd_SE_cat['xpeak_image'] = XPEAK_IMAGE
			ccd_SE_cat['ypeak_image'] = YPEAK_IMAGE
			ccd_SE_cat['x_image'] = X_IMAGE
			ccd_SE_cat['y_image'] = Y_IMAGE
			ccd_SE_cat['RA'] = ALPHA_J2000
			ccd_SE_cat['DEC'] = DELTA_J2000	
			
			DWF_X = Table()
			DWF_X['RA'] = ALPHA_J2000
			DWF_X['DEC'] = DELTA_J2000
			
			RA, DEC, Gmag_vega, Gmag_vega_err, Gmag_AB = np.loadtxt(output_cats + 'gaiaDR2.ascii', unpack = True)
			print(RA) 



			'''
			#Arecno, ARAJ2000, ADEJ2000, Ae_RAJ2000, Ae_DEJ2000 ,AField, Anobs, Amobs, AB_V, Ae_B-V, AVmag, Ae_Vmag, ABmag, Ae_Bmag, Ag_mag, Ae_g_mag, Ar_mag, Ae_r_mag, Ai_mag, Ae_i_mag = np.loadtxt(output_cats + 'AAVOS.cat', unpack=True)
			Arecno, ARAJ2000, ADEJ2000, Ae_RAJ2000, Ae_DEJ2000, AField, Anobs, Amobs, AB_V, Ae_B, AVmag, Ae_Vmag, ABmag, Ae_Bmag, Ag_mag, Ae_g_mag, Ar_mag, Ae_r_mag, Ai_mag, Ae_i_mag = np.loadtxt(output_cats + 'AAVOS.cat', unpack=True)
			print(Arecno)
			AAVOS_cat = Table()
			AAVOS_cat['RA'] = ARAJ2000
			AAVOS_cat['DEC'] = ADEJ2000
			AAVOS_cat['g_mag'] = Ag_mag
			AAVOS_cat['g_mag_err'] = Ae_g_mag
			AAVOS_cat['r_mag'] = Ar_mag
			AAVOS_cat['r_mag_err'] = Ae_r_mag
			AAVOS_cat['i_mag'] = Ai_mag
			AAVOS_cat['i_mag_err'] = Ae_i_mag
			
			AAVOS_X = Table()
			AAVOS_X['RA']= AAVOS_cat['RA'] 
			AAVOS_X['DEC']= AAVOS_cat['DEC']
			print(AAVOS)
			
			#USNO-B1, URAJ2000, UDEJ2000, Ue_RAJ2000, Ue_DEJ2000, UEpoch, UpmRA, UpmDE ,UNdet ,UB1mag, UR1mag ,UB2mag ,UR2mag, UImag = np.loadtxt(output_cats + USNO_B1_output, unpack=True)
			USNO, URAJ2000, UDEJ2000, Ue_RAJ2000, Ue_DEJ2000, UEpoch, UpmRA, UpmDE, UNdet, UB1mag, UR1mag, UB2mag, UR2mag, UImag = np.loadtxt(output_cats + 'USNO_B1.cat', unpack=True)
			USNO_B1_cat = Table()
			USNO_B1_cat['RA'] = URAJ2000
			USNO_B1_cat['DEC'] = UDEJ2000
			USNO_B1_cat['R1_mag'] = UR1mag
			USNO_B1_cat['R2_mag'] = UR2mag
			USNO_B1_cat['I_mag'] = UImag
			
			USNO_B1_X = Table()
			USNO_B1_X['RA'] = URAJ2000
			USNO_B1_X['DEC'] = UDEJ2000
			'''
			
			##---Cross match with gaiaDR2 ------ # 
			'''
			max_radius = 1./3600 #1 arc second 
			
			dist_between, ind_row = crossmatch_angular(DWF_X, AAVOS_X, max_radius)
			match = ~np.isinf(dist_between)
			
			match_table_aavos = Table()
			match_table_aavos['matched_true_false'] = match 
			match_table_aavos['matched_ID'] = ind_row
			match_table_aavos['matched_DWF_data_gmag'] = AAVOS_cat['g_mag']
			match_table_aavos['matched_DWF_data_gmag_err'] = AAVOS_cat['g_mag_err']
			match_table_aavos['matched_DWF_RA'] = AAVOS_cat['RA']
			match_table_aavos['matched_DWF_DEC'] = AAVOS_cat['DEC']
			
			print(match_table_aavos)	
			'''
