from astropy.table import Table
from svnds import PATH_DATA
import os

# import LSST filters
path_filters = os.path.join(PATH_DATA, 'filters/LSST/')
u_L = Table.read(path_filters + 'total_u.dat', format = 'ascii')
g_L = Table.read(path_filters + 'total_g.dat', format = 'ascii')
r_L = Table.read(path_filters + 'total_r.dat', format = 'ascii')
i_L = Table.read(path_filters + 'total_i.dat', format = 'ascii')
z_L = Table.read(path_filters + 'total_z.dat', format = 'ascii')
y_L = Table.read(path_filters + 'total_y.dat', format = 'ascii')


# Parameters
FILTERS = [u_L, g_L, r_L, i_L, z_L, y_L]
NFILTER = len(FILTERS)
NAMES = ['u_L','g_L', 'r_L', 'i_L', 'z_L', 'y_L']
MAG5 = [26.1, 27.4, 27.5, 26.8, 26.1, 24.9] # S/R = 5 limiting magnitudes from Graham 2020 and Ivezi 2019

COL_WAVE = 'col1'
FACTOR_WAVE = 1 / 1000 # nm to microns
COL_RESP = 'col2'