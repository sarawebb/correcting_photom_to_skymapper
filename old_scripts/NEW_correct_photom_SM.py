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
input_path = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/*/SE_cats/'
output_path = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/'
input_cats = '/mnt/dwf/sky_mapperDR2/shortlisted_stars_for_photom/4hr_SM.csv_SHORTLISTED.ascii'

###-----------------------------------^^^^^^^^^^^ DID YOU INPUT FIELD NAMEs ? ^^^^^^ --------------------------------------------------------------


path_list = glob.glob(input_path)

for i in path_list: 
	for filename in os.listdir(i):
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
			
			#print(ccd_SE_cat)
			
			os.chdir(i)
			os.chdir("..")
			print(os.path.abspath(os.curdir))
			
			#photom_correction_path = 
			
