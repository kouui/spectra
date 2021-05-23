


from spectra_src.ImportAll import *

from spectra_src.Struct import Atom, Atmosphere, Radiation
from spectra_src.Function.SEquil import SELib
from spectra_src.Function import SlabModel 

import os

conf_path = os.path.join( CFG._ROOT_DIR, "data/conf/H.conf" )
atom, wMesh, path_dict = Atom.init_Atom_(conf_path,is_hydrogen=True)

atmos = Atmosphere.Atmosphere0D(Nh=1.E12, Ne=1.E11, Te=7.E3, Vd=0., Vt=5.E5)
radiation = Radiation.init_Radiation_(atmos, wMesh, 0.5)
SE_con, Rate_con = SELib.cal_SE_with_Nh_Te_(atom, atmos, wMesh,radiation, Nh_SE = None)
Cloud_con = SlabModel.SE_to_slab_0D_(atom, atmos, SE_con, depth = 1.E3 * 1.E5) # 1_000 [km]

# Cloud_con.w0[:]             central wavelength in [cm]
# Cloud_con.tau_max[:]        maximum optical depth
# Cloud_con.Ibar[:]           integration of (intensity_profile * wavelength_mesh)

print(f"{'WAVELENGTH [AA]'}  {'MAX TAU'}      {'INTEGRATED INTENSITY [erg/cm^2/Sr/s]'}")
for k in (7,8,9):
    print(f"{Cloud_con.w0[k]*1.E8:.2f}{' ':10s}{Cloud_con.tau_max[k]:.2E}{' ':5s}{Cloud_con.Ibar[k]:.2E}")
