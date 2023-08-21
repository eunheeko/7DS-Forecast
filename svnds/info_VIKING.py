from astropy.table import Table

# import EUCLID filters
path_filters = './data/filters/VIKING/'
Z_V = Table.read(path_filters + 'VISTA_Filters_at80K_forETC_Z.dat', format = 'ascii.no_header')
Y_V = Table.read(path_filters + 'VISTA_Filters_at80K_forETC_Y.dat', format = 'ascii.no_header')
J_V = Table.read(path_filters + 'VISTA_Filters_at80K_forETC_J.dat', format = 'ascii.no_header')
H_V = Table.read(path_filters + 'VISTA_Filters_at80K_forETC_H.dat', format = 'ascii.no_header')
Ks_V = Table.read(path_filters + 'VISTA_Filters_at80K_forETC_Ks.dat', format = 'ascii.no_header')



# Parameters
FILTERS = [Z_V, Y_V, J_V, H_V, Ks_V]
NFILTER = len(FILTERS)
NAMES = ['Z_V', 'Y_V', 'J_V', 'H_V', 'Ks_V']
MAG5 = [23.1, 22.3, 22.1, 21.5, 21.2] # S/R = 5 limiting magnitudes from VKING DR4

COL_WAVE = 'col1'
FACTOR_WAVE = 1 / 1000 # nm to micron
COL_RESP = 'col2'