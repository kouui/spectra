
#-------------------------------------------------------------------------------
# uiltility to manipulate VERSION related issue in spectra
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

import os, sys
import re
from typing import List


def scan_version_(dir : str) -> None:
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
                get_file_version_(f"{cur_dir}/{file}")

    return None


def get_file_version_( fname : str ) -> None:
    """print out {VERSION} and {fname} in a single lline

    Args:
        fname (str): full path to the file

    Returns:
        None: None
    """
    is_version = False
    is_found   = False
    with open(fname, 'r') as f:
        for line in f:

            if line.startswith("# VERSION"):
                is_version = True
                continue
            
            if is_version:
                words : List[str] = [w.strip() for w in line.split()]
                for word in words:
                    match = re.match(r'[0-9].[0-9].[0-9]', word)
                    
                    if match is None:
                        continue
                    
                    print(f"{match.string:<8s}   {fname}")
                    is_found = True

            if is_found:
                break
        
        if not is_found:
            print(f"{'None':<8s}   {fname}")
    return None


if __name__ == "__main__":

    print('-'*50)
    print(f"{'VERSION':<8s}   {' '*1}{'FILE'}")
    print('-'*50)

    dir : str = sys.argv[1]
    scan_version_( dir )