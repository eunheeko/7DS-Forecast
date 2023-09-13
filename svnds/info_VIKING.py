from astropy.table import Table
from svnds import PATH_DATA
import os


# import VIKING filters
path_filters = os.path.join(PATH_DATA, 'filters/VIKING/')
Z_V = Table.read(path_filters + 'VISTA_Filters_at80K_forETC_Z.dat', format = 'ascii.no_header')
Y_V = Table.read(path_filters + 'VISTA_Filters_at80K_forETC_Y.dat', format = 'ascii.no_header')
J_V = Table.read(path_filters + 'VISTA_Filters_at80K_forETC_J.dat', format = 'ascii.no_header')
H_V = Table.read(path_filters + 'VISTA_Filters_at80K_forETC_H.dat', format = 'ascii.no_header')
Ks_V = Table.read(path_filters + 'VISTA_Filters_at80K_forETC_Ks.dat', format = 'ascii.no_header')



# Parameters
FILTERS = [Z_V, Y_V, J_V, H_V, Ks_V]
NFILTER = len(FILTERS)
NAMES = ['Z_V', 'Y_V', 'J_V', 'H_V', 'Ks_V']

'''
From Barnett+21, 5-sigma limiting magnitudes at (ZYJHKs) is (22.1, 21.3, 20.9, 19.8, 19.2) in VEGA
    - NOTE: It seems that the paper has typo (22.1, `20.3`, 20.9, 19.8, 19.2) but MAG5 at Y mag is 21.3 according to etheFigure 1

However, from VIKING DR4 paper in Venemans+13 (Table1), (23.1, 22.3, 22.1, 21.5, 21.2) in AB mag
- This uses conversion eq as (0.521, 0.618, 0.937, 1.384, 1.839)

Here, we adopt the value from Barnett+21 and conversion equation from (http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/filter-set)
for two reasons:
    (1) Barnett+21 is more up-to-date
    (2) Total survey area covers a larger region in Barnett

Vega to AB conversions are (0.502, 0.600, 0.916, 1.366, 1.827) (v1.3 ver)

'''
MAG5 = [22.1 + 0.502, 
        21.3 + 0.600, 
        20.9 + 0.916, 
        19.8 + 1.366,
        19.2 + 1.827
        ]
# MAG5 = [23.1, 22.3, 22.1, 21.5, 21.2] # S/R = 5 limiting magnitudes from VKING DR4


COL_WAVE = 'col1'
FACTOR_WAVE = 1 / 1000 # nm to micron
COL_RESP = 'col2'