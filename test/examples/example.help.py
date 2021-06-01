

from spectra_src.ImportAll import *

from spectra_src.Struct import Atom
from spectra_src.Util import HelpUtil

import os

conf_path = os.path.join( CFG._ROOT_DIR, "data/conf/H.conf" )
atom, wMesh, path_dict = Atom.init_Atom_(conf_path,is_hydrogen=True)

HelpUtil.help_( atom )