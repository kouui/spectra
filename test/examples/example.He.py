

from spectra_src.ImportAll import *

from spectra_src.Struct import Atom, Atmosphere, Radiation
from spectra_src.Function.SEquil import SELib

import os

conf_path = os.path.join( CFG._ROOT_DIR, "data","conf","He.conf" )
atom, wMesh, path_dict = Atom.init_Atom_(conf_path,is_hydrogen=False)
breakpoint()
atmos = Atmosphere.Atmosphere0D(Nh=1.E12, Ne=1.E11, Te=7.E3, Vd=0., Vt=5.E5)
radiation = Radiation.init_Radiation_(atmos, wMesh, 0.5)
SE_con, Rate_con = SELib.cal_SE_with_Nh_Te_(atom, atmos, wMesh,radiation, None)

print(f">>> calculate SE given Nh and Te <<<\n")
print(f"Electron temperature  = {atmos.Te:.1E}")
print(f"Electron density      = {atmos.Ne:.1E}")
print('-'*35)
print("SE  :")
for v in SE_con.n_SE[:]:
    print(f"{v:.4E}",end="  ")
print(f"\nLTE :")
for v in SE_con.n_LTE[:]:
    print(f"{v:.4E}",end="  ")
print(" ")


"""Correct answer
"""

atmos = Atmosphere.Atmosphere0D(Nh=1.E11, Ne=5.E10, Te=7.E3, Vd=0., Vt=5.E5)
radiation = Radiation.init_Radiation_(atmos, wMesh, 0.5)
SE_con, Rate_con = SELib.cal_SE_with_Ne_Te_(atom, atmos, wMesh,radiation, None)

print(f"\n>>> calculate SE given Ne and Te <<<\n")
print(f"Electron temperature  = {atmos.Te:.1E}")
print(f"Hydrogen density      = {atmos.Nh:.1E}")
print('-'*35)
print("SE  :")
for v in SE_con.n_SE[:]:
    print(f"{v:.4E}",end="  ")
print(f"\nLTE :")
for v in SE_con.n_LTE[:]:
    print(f"{v:.4E}",end="  ")
print(" ")