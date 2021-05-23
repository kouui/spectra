

from spectra_src.ImportAll import *

from spectra_src.Struct import Atom, Atmosphere, Radiation
from spectra_src.Function.SEquil import SELib

import os

conf_path = os.path.join( CFG._ROOT_DIR, "data/conf/H.conf" )
atom, wMesh, path_dict = Atom.init_Atom_(conf_path,is_hydrogen=True)

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
Electron temperature  = 7.0E+03
Electron density      = 1.3E+11
-----------------------------------
SE  :
8.6628E-01  3.3842E-07  3.0978E-09  1.1299E-09  8.7587E-10  8.3796E-10  8.9827E-10  1.0387E-09  1.3372E-01  
LTE :
3.7019E-01  6.7241E-08  6.6075E-09  3.9264E-09  3.6943E-09  4.0386E-09  4.6556E-09  5.4592E-09  6.2981E-01 
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