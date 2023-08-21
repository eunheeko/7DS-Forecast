from svnds.synphot import Syn7DS, merge_tbl
import numpy as np
import os

mag_ref = 23
path_elcosmos = '/data8/EL_COSMOS/elcosmos/'
path_save = f'/data8/EL_COSMOS/synphots/SDS/scale_{mag_ref:d}_eff/'

Tsamps = np.array([180, 180 * 364/14, 180 * 364/14 * 3, 180 * 364/14 * 5, 180 * 364/14 * 7, 180 * 365, 180 * 365 * 5], dtype = int)
Tsamps = (Tsamps * 0.5).astype(int)

syn = Syn7DS(path_elcosmos, path_save, Tsamps)

syn.synphot(zones = range(1, 10), magscale = True, mag_ref = mag_ref)

merge_tbl(path_save, f"synphot_7DS_scale_{mag_ref:d}_yrs_all_eff_tot.csv")

