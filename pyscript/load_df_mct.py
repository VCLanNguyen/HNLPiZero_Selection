import sys
import os

os.nice(20)

import warnings
warnings.filterwarnings("ignore")

import argparse
import uproot
import h5py
import pandas as pd
import time

# Local helper script
hnlDIR = os.environ['_']
sys.path.append(hnlDIR + '/pyscript/')

from Branches import *
from HelperFunctions import *
from CutFunctions import *

def load_tree(path, tname, branches):
    
    folder = "hnlana"
    tree = uproot.open(path = "{}:{}/{}".format(path, folder, tname))
                        #, object_cache=8000
                        #, array_cache = "8000 MB"
                       #, xrootd_handler=uproot.MultithreadedXRootDSource
                       #, num_workers=10)
			
    df = tree.arrays(branches, library = "pd")

    return df

def load_df(input_file):
        
    dfmct = load_tree(input_file, "events", mct_branches)

    return dfmct

def make_new_dataframe():
    
    dfmct_cc  = pd.DataFrame() 

    return dfmct_cc

def reset_index_dataframe(dfmct_cc):
    dfmct_cc = dfmct_cc.set_index(['run','subrun', 'event']).sort_index().reset_index()

    return dfmct_cc

def main(args):

    file_list = open(args.input, "r")
    lines = file_list.readlines()
    file_list.close()
    nlines = len(lines) - 1
    
    start = time.time()
    
    dfmct_cc = make_new_dataframe()
    itr = 0
    
    for i, line in enumerate(lines):

        print("Reading in files " + str(i) + ": " + str(line))
        dfmct = load_df(line)
   
        dfmct_cc = pd.concat([dfmct_cc, dfmct], ignore_index = True) 

        if ((i != 0) and (i%40 == 0)) or ((i != 0) and (i%nlines == 0)) :
    
            print ("Saving iteration " + str(itr) + ", file counter " + str(i))
            
            reset_index_dataframe(dfmct_cc)       
    
            ##save to pickle
            dfmct_cc.to_pickle(args.output+"_mct_" + str(itr) + ".pkl", protocol = 5)
    
            del dfmct_cc
            
            ##make new dataframe and restart the cycle
            dfmct_cc = make_new_dataframe()
            itr += 1

    now = time.time()
    duration = (now - start)/3600
    print("Done saving root to dataframe, elapsed time = " + str(duration) + " hour(s)") 

if __name__ == '__main__':

    parser = argparse.ArgumentParser()      
    parser.add_argument("--input", required=True, help="Input list of root file e.g. hnl.list")
    parser.add_argument("--output", required=True, help="Output path for pkl file e.g. /save/to/this/path/")
    parser.add_argument("--type", required=True, help="hnl or nu or cos?")

    args = parser.parse_args()                             
    main(args)                                             
