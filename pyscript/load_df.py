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
    
    folder = "hnlpizeroana"
    tree = uproot.open(path = "{}:{}/{}".format(path, folder, tname)
                        , object_cache=2000
                        , array_cache = "2000 MB")
    df = tree.arrays(branches, library = "pd")

    return df

def load_df(input_file, output_file, ifSlc = True, ifSubrun = True, ifMeVPrtl = False, ifMct = False ):
        
    #SLICE AND PFP NEED SPECIAL HANDLING
    dfslc = load_tree(input_file, "events", slc_branches)
    dfpfp = load_tree(input_file, "events", pfp_branches)
   
    #add slc idx for joining
    dfslc["slc_idx"] = dfslc.groupby(['run', 'subrun', 'event']).transform("cumcount").add(1) - 1
    #every row = slice
    dfslc = dfslc.set_index(['run', 'subrun', 'event']).reset_index()
    
    #explode every row = slice
    dfpfp = dfpfp.set_index(['run', 'subrun', 'event']).apply(pd.Series.explode).reset_index()
    #change stlvector to array
    dfpfp = dfpfp.set_index(['run', 'subrun', 'event']).applymap(lambda x: np.array(x))
    #add slc idx for joining
    dfpfp["slc_idx"] = dfpfp.groupby(['run', 'subrun', 'event']).transform("cumcount").add(1) - 1
    dfpfp = dfpfp.reset_index()
   
    #merge pfp and slc dataframe on slice idx
    dfjoin = pd.concat([dfslc, dfpfp], axis = 1)
    dfjoin = dfjoin.loc[:, ~dfjoin.columns.duplicated()]
    #explode to pfp level
    dfjoin = dfjoin.apply(pd.Series.explode) #explode to pfp level
    
    #reset index? 
    #dfjoin = dfjoin.set_index(['run', 'subrun', 'event', 'slc_idx']).sort_index().reset_index()
    #dfjoin = dfjoin.set_index(['run', 'subrun', 'event', 'slc_idx']).sort_index()

    #Pre-selection cut
    #dfjoin = cutPreSelection(dfjoin)

    #CRUMBS cut
    #dfjoin = cutCosmics(dfjoin)
    
    #fix unit
    dfjoin["slc_opt0_time_corrected_Z_pandora"] = dfjoin["slc_opt0_time_corrected_Z_pandora"] * 1000
    dfjoin.to_hdf(output_file, key='slc', mode = 'w')
    
    dfsubrun = load_tree(input_file, "subruns", subrun_branches)
    dfsubrun = dfsubrun.set_index(['run','subrun']).sort_index().reset_index()
    dfsubrun.to_hdf(output_file, key='subrun')
    
    dfmevprtl = load_tree(input_file, "events", mevprtl_branches)
    dfmevprtl = dfmevprtl.set_index(['run','subrun', 'event']).sort_index().reset_index()
    dfmevprtl.to_hdf(output_file, key='mevprtl')
    
    dfmct = load_tree(input_file, "events", mct_branches)
    dfmct = dfmct.set_index(['run','subrun', 'event']).sort_index().reset_index()
    dfmct.to_hdf(output_file, key='mct')

def main():

    hnl_input_file = "./root_files/hnl_m200_20k.root"
    hnl_output_file = "./hdf5_files/hnl_bdt_test_20k.h5"

    nu_input_file = "./root_files/tpc_nu_20k.root"
    nu_output_file = "./hdf5_files/nu_bdt_test_20k.h5"
    
    cosmics_input_file = "./root_files/cosmics_5k.root"
    cosmics_output_file = "./hdf5_files/cosmics_test_5k.h5"
    
    start = time.time()

    load_df(hnl_input_file, hnl_output_file, ifSlc = True, ifSubrun = True, ifMeVPrtl = False, ifMct = False)
    
    now1 = time.time()
    print("Done saving hnl, elapsed time = " + str(now1 -start))
    
    load_df(nu_input_file, nu_output_file, ifSlc = True, ifSubrun = True, ifMeVPrtl = False, ifMct = False)
    
    now2 = time.time()
    print("Done saving nu, elapsed time = " + str(now2 - now1))
    
    #load_df(cosmics_input_file, cosmics_output_file, ifSlc = True, ifSubrun = True, ifMeVPrtl = False, ifMct = False)

    #now3 = time.time()
    #print("Done saving cosmics, elapsed time = " + str(now3 - now2))

if __name__ == '__main__':

	main()
