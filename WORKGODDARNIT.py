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

##-----------------------------------INPUT FIELD NAME --------------------------------------------------------------


input_path = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/c4d_150116_053736_ooi_g_v1/SE_cats/'
output_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'
input_cats = '/mnt/dwf/sky_mapperDR2/4hr_SM.csv'

path_list = glob.glob(input_path)

for i in path_list:
    for filename in os.listdir(i):
        if filename.endswith('c4d_150116_053736_ooi_g_v1._ext_10.fits_thresh_1.5_SE.cat'):
        
            MAG_APER, MAGERR_APER, MAG_AUTO, MAGERR_AUTO, XPEAK_IMAGE, YPEAK_IMAGE, X_IMAGE,Y_IMAGE, ALPHA_J2000, DELTA_J2000 = np.loadtxt(i + filename, unpack = True)

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

            DWF_X = np.empty((len(ccd_SE_cat), 2), dtype=np.float64)
            DWF_X[:, 0] = ALPHA_J2000
            DWF_X[:, 1] = DELTA_J2000

            SM_cat  = pd.read_csv(input_cats).to_dict(orient='row')
        
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
                if 0.80 <= row['class_star'] <= 1.0:
                    if 4 <= row['ngood'] <= 20:
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

            SM_X = np.empty((len(SM_filtered), 6), dtype=np.float64 )
            SM_X[:,0] = SM_filtered['ra']
            SM_X[:,1] = SM_filtered['dec']
            SM_X[:,2] = SM_filtered['g_psf']
            SM_X[:,3] = SM_filtered['r_psf']
            SM_X[:,4] = SM_filtered['i_psf']
            SM_X[:,5] = SM_filtered['z_psf']
            print(SM_X)
            ##---------------------------------Cross match with SkyMapper -------#
			max_radius = 2./3600 #2 arc seconds 
            dist_between, ind_row = crossmatch_angular(DWF_X, SM_X, max_radius)
            match = ~np.isinf(dist_between)
            
			match_table = Table()
            match_table['matched_true_false'] = match 
            match_table['matched_ID'] = ind_row
            match_table['matched_DWF_data_gmag'] = ccd_SE_cat['mag_auto'] = MAG_AUTO
            match_table['matched_DWF_RA'] = ccd_SE_cat['RA']
            match_table['matched_DWF_DEC'] = ccd_SE_cat['DEC']
            print(match_table)
			
            SM_match_true = []
            SM_row_matched = []
            DWF_g_mags_matched = []
            DWF_g_mags_error_matched = []
            DWF_obs_ra_matched = []
            DWF_obs_dec_matched = []
            
            for row in match_table: 
                if row['matched_true_false'] == True: 
                    SM_match_true.append(row['match_true_false'])
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
            for i in SM_match_true:
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
            print(source_match_table)
			
			
