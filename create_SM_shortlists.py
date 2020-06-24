
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
from astroML.crossmatch import crossmatch_angular 
from astroML.datasets import fetch_imaging_sample, fetch_sdss_S82standards
from astroML.plotting import hist


input_cats = '/fred/oz100/NOAO_archive/archive_NOAO_data/scripts/correct_photom/skymapper/'

cat_list = glob.glob(input_cats)

for filename in os.listdir(input_cats):
		if filename.endswith('.csv'):

			SM_cat  = pd.read_csv(input_cats + filename).to_dict(orient='row')
			#print(SM_cat[0]['i_psf'])



			SM_ra = []
			SM_dec = []
			SM_ngood = []
			SM_class_star = []
			SM_g_psf = []
			SM_g_psf_err =[]
			SM_r_psf = []
			SM_r_psf_err =[]
			SM_i_psf = []
			SM_i_psf_err =[]
			SM_z_psf = []
			SM_z_psf_err =[]

			for i in range(len(SM_cat)):
    				SM_ra.append(SM_cat[i]['raj2000'])
    				SM_dec.append(SM_cat[i]['dej2000'])
    				SM_ngood.append(SM_cat[i]['ngood'])
    				SM_class_star.append(SM_cat[i]['class_star'])
    				SM_g_psf.append(SM_cat[i]['g_psf'])
    				SM_g_psf_err.append(SM_cat[i]['e_g_psf'])
    				SM_r_psf.append(SM_cat[i]['r_psf'])
    				SM_r_psf_err.append(SM_cat[i]['e_r_psf'])
    				SM_i_psf.append(SM_cat[i]['i_psf'])
    				SM_i_psf_err.append(SM_cat[i]['e_i_psf'])
    				SM_z_psf.append(SM_cat[i]['z_psf'])
    				SM_z_psf_err.append(SM_cat[i]['e_z_psf'])
    
			SM_table = Table()
			SM_table['RA'] = SM_ra
			SM_table['DEC'] = SM_dec
			SM_table['ngood'] = SM_ngood
			SM_table['class_star'] = SM_class_star
			SM_table['g_psf'] = SM_g_psf
			SM_table['e_g_psf'] = SM_g_psf_err
			SM_table['r_psf'] = SM_r_psf
			SM_table['e_r_psf'] = SM_r_psf_err
			SM_table['i_psf'] = SM_i_psf
			SM_table['e_i_psf'] = SM_i_psf_err
			SM_table['z_psf'] = SM_z_psf
			SM_table['e_z_psf'] =SM_z_psf_err

			#print(SM_table)


			SM_filered_ra = []
			SM_filered_dec = []
			SM_filered_ngood = []
			SM_filered_class_star = []
			SM_filered_g_psf = []
			SM_filered_g_psf_err = []
			SM_filered_r_psf = []



			SM_filered_r_psf_err = []
			SM_filered_i_psf = []
			SM_filered_i_psf_err = []
			SM_filered_z_psf = []
			SM_filered_z_psf_err = []
			
			for row in SM_table:
				if 0.95 <= row['class_star'] <= 1.0:
					SM_filered_ra.append(row['RA'])
					SM_filered_dec.append(row['DEC'])
					SM_filered_ngood.append(row['ngood'])
					SM_filered_class_star.append(row['class_star'])
					SM_filered_g_psf.append(row['g_psf'])
					SM_filered_g_psf_err.append(row['e_g_psf'])
					SM_filered_r_psf.append(row['r_psf'])
					SM_filered_r_psf_err.append(row['e_r_psf'])
					SM_filered_i_psf.append(row['i_psf'])
					SM_filered_i_psf_err.append(row['e_i_psf'])
					SM_filered_z_psf.append(row['z_psf'])
					SM_filered_z_psf_err.append(row['e_z_psf'])
			
			
			
			'''if row in SM_table is 0.80 <= row['class_star'] <= 1.0 and 5 <= row['ngood'] <= 20:
				SM_filered_ra.append(row['RA'])
				SM_filered_dec.append(row['DEC'])
				SM_filered_ngood.append(row['ngood'])
				SM_filered_g_psf.append(row['g_psf'])
				SM_filered_g_psf_err.append(row['e_g_psf'])
				SM_filered_r_psf.append(row['r_psf'])
				SM_filered_r_psf_err.append(row['e_r_psf'])
				SM_filered_i_psf.append(row['i_psf'])
				SM_filered_i_psf_err.append(row['e_i_psf'])
				SM_filered_z_psf.append(row['z_psf'])
				SM_filered_z_psf_err.append(row['e_z_psf'])
			
			for row in SM_table: 
    				if 0.80 <= row['class_star'] <= 1.0 and 5 <= row['ngood'] <= 20: 
                        SM_filered_ra.append(row['RA'])
                        SM_filered_dec.append(row['DEC'])
                        SM_filered_ngood.append(row['ngood'])
                        SM_filered_class_star.append(row['class_star'])
                        SM_filered_g_psf.append(row['g_psf'])
                        SM_filered_g_psf_err.append(row['e_g_psf'])
                        SM_filered_r_psf.append(row['r_psf'])
                        SM_filered_r_psf_err.append(row['e_r_psf'])
                        SM_filered_i_psf.append(row['i_psf'])
                        SM_filered_i_psf_err.append(row['e_i_psf'])
                        SM_filered_z_psf.append(row['z_psf'])
                        SM_filered_z_psf_err.append(row['e_z_psf'])
	    		'''	
				
			SM_filtered = Table()
			SM_filtered['ra'] =  SM_filered_ra 
			SM_filtered['dec'] =  SM_filered_dec
			SM_filtered['ngood'] =  SM_filered_ngood
			SM_filtered['class_star'] =  SM_filered_class_star
			SM_filtered['g_psf'] =  SM_filered_g_psf 
			SM_filtered['e_g_psf'] =  SM_filered_g_psf_err 
			SM_filtered['r_psf'] =  SM_filered_r_psf 
			SM_filtered['e_r_psf'] =  SM_filered_r_psf_err 
			SM_filtered['i_psf'] =  SM_filered_i_psf 
			SM_filtered['e_i_psf'] =  SM_filered_i_psf_err 
			SM_filtered['z_psf'] =  SM_filered_z_psf 
			SM_filtered['e_z_psf'] =  SM_filered_z_psf_err 
			
			t = SM_filtered
			print(t)
			output = '/fred/oz100/NOAO_archive/archive_NOAO_data/scripts/correct_photom/skymapper/shortlisted_stars_for_photom/' + filename + '_SHORTLISTED.ascii'
			t.write(output, format='ascii', overwrite= True)
			
