# mypy: ignore-errors

#-------------------------------------------------------------------------------
# configurations in spectra
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------


from .Types import T_STR, T_BOOL, T_INT, T_NORETURN
import os

_ROOT_DIR : T_STR = os.path.abspath( os.path.join( os.path.dirname(__file__), ".." ) )


#-------------------------------------------------------------------------------
# for numba configuratoin
#-------------------------------------------------------------------------------



_IS_JIT       : T_BOOL = False

_IS_NOGIL     : T_BOOL = False

_IS_PARALLEL  : T_BOOL = False

_IS_CACHE     : T_BOOL = False

"""whether to turn on the JIT compilation in all *.py files,
since sphinx does no understand numba thus can not generate documentation for numba jitted functions.

set to
    True  : before simulation; before pushing to github

    False : before generating documentation using sphinx
"""

from numba.core import config as nb_config
from numba import set_num_threads as nb_set_num_threads

def _SET_NUMBA_THREAD_(threading_layer : T_STR ='threadsafe', n_thread : T_INT = 2) -> None:
    r""" """

    nb_config.THREADING_LAYER = threading_layer
    """ set to thread safe library : tbb, omp, workqueue, default : workqueue. tbb if possible """

    nb_set_num_threads(n_thread)
    """ limiting the number of threads """

    return None

