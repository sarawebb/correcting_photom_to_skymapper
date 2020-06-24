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
input_path = '/mnt/dwf/archive_NOAO_data/data_outputs/*/*/*/g_band/single/*/photom_correction_files/skymapper_checkphotom_plots/'


###-----------------------------------^^^^^^^^^^^ DID YOU INPUT FIELD NAMEs ?  --------------------------------------------------------------

path_list = glob.glob(input_path)
#print(path_list)

for i in path_list: 
	print(i)
	for filename in os.listdir(i):
		#av_correction = np.zeros(len(os.listdir(i)))
		#print(av_correction)
		filenames_for_av = []
		filenames_for_av.append(filename)
		if filename.endswith('.png'):
			os.system('rm ' + str(i) + str(filename))
			print('removing ' + str(filename)) 
				
			
		 
			
