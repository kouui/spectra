

#-------------------------------------------------------------------------------
# check the type of all *.py file in the directory recursively using mypy 
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

import os
import sys
#from mypy.__main__ import console_entry

_FILEs_IGNORE = {
    "Configurations.py" : None
}

def typecheck_folder_(dir : str) -> None:
    """run get_file_version_ for all selected files in the dir, recursively

    Args:
        dir (str): path to the directory

    Returns:
        None: None
    """

    for cur_dir, in_dir, files in os.walk(dir):
        
        ## : skip ".*" folder
        if os.path.basename( cur_dir ).startswith('.'):
            continue

        for file in files:
            ## : skip "_*" file
            if file.startswith('_'):
                continue
            ## : skip "Removal.*" file
            if file.startswith("Removal."):
                continue
            ## : pick "*.py" file
            if file.endswith(".py"):
                print( f"{cur_dir}/{file}", end='' )
                try : 
                    _ = _FILEs_IGNORE[file]
                    print(" --> ignore")
                    continue
                except KeyError:
                    print('')
                
                
                os.system( f"mypy {cur_dir}/{file}" )
                #sys.argv[0] = f"{cur_dir}/{file}"
                #console_entry()

    return None

if __name__ == "__main__":

    typecheck_folder_(sys.argv[1])

