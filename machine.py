from svnds.synphot import Syn7DS, SynSPx, merge_tbl
import numpy as np
import os
from astropy.table import Table

for mag_ref in [20.5, 21.5, 22.5]:

    path_elcosmos = '/data8/EL_COSMOS/elcosmos/'
    # path_save = f'/data8/EL_COSMOS/synphots/SDS_m6250/scale_{mag_ref:.1f}_eff/'
    path_save = f'/data8/EL_COSMOS/synphots/SPHEREx_m6250/scale_{mag_ref:.1f}/'

    # Tsamps = np.array([180 * 364/14 * 5])
    # Tsamps = (Tsamps * 0.5).astype(int)
    syn = SynSPx(path_elcosmos, path_save)

    syn.synphot(zones = range(1, 10), magscale = True, mag_ref = mag_ref)

    merge_tbl(path_save, f"synphot_SPHEREx_scale_{mag_ref:.1f}_tot.csv")
