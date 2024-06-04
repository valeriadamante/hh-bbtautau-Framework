import ROOT
import sys
import os
import math
import shutil
from RunKit.run_tools import ps_call
if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])

import Common.Utilities as Utilities
from Analysis.HistHelper import *

if __name__ == "__main__":
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFile', nargs='+', type=str)
    parser.add_argument('--outFile', required=True)
    parser.add_argument('--var', required=False, type=str, default='tau1_pt')
    #parser.add_argument('--remove-files', required=False, type=bool, default=False)
    args = parser.parse_args()

    # 1 list files :


    all_files = [ fileName for fileName in args.inputFile ]
    hadd_str = f'hadd -f209 -n 10 {args.outFile} '
    hadd_str += ' '.join(f for f in all_files)
    print(hadd_str)
    if len(all_files) > 1:
        ps_call([hadd_str], True)
    else:
        shutil.copy(all_files[0],args.outFile)
    #if os.path.exists(outFileNameFinal) and args.remove_files:
    #    for histFile in all_files[var]:
    #        if histFile == outFileNameFinal: continue
    #        os.remove(histFile)