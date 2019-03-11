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

input_path = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/Prime/g_band/single/*/photom_correction_files/'

path_list = glob.glob(input_path)

for i in path_list: 
	averages = []
	for filename in os.listdir(i):
		if filename.startswith('average'):
			try: 
				avs = np.loadtxt(i + filename, unpack = True, skiprows= 1)
				#print(avs)
				for f in avs: 
					averages.append(f)
			except:
				pass 
	#print(averages)
	print(i)
	print(np.mean(averages))
