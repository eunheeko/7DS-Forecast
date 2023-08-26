from astropy.table import Table
from svnds import PATH_DATA
import os

# import EUCLID filters
path_filters = os.path.join(PATH_DATA, 'filters/EUCLID/')
J_E = Table.read(path_filters + 'NISP-PHOTO-PASSBANDS-V1-J_throughput.dat', format = 'ascii.commented_header')
H_E = Table.read(path_filters + 'NISP-PHOTO-PASSBANDS-V1-H_throughput.dat', format = 'ascii.commented_header')
Y_E = Table.read(path_filters + 'NISP-PHOTO-PASSBANDS-V1-Y_throughput.dat', format = 'ascii.commented_header')



# Parameters
FILTERS = [Y_E, J_E, H_E]
NFILTER = len(FILTERS)
NAMES = ['Y_E', 'J_E', 'H_E']
MAG5 = [24.0, 24.2, 23.9] # S/R = 5 limiting magnitudes  from Graham et al. 2020, Tssable1

COL_WAVE = 'WAVE'
FACTOR_WAVE = 1 / 1000 # nm to micron
COL_RESP = 'T_TOTAL'