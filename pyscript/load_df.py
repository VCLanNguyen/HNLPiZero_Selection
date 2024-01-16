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

def load_df(input_file, output_file, ftype):
        
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

    del dfpfp
    del dfslc

    dfjoin = dfjoin.loc[:, ~dfjoin.columns.duplicated()]
    #explode to pfp level
    dfjoin = dfjoin.apply(pd.Series.explode) #explode to pfp level
    
    #reset index? 
    #dfjoin = dfjoin.set_index(['run', 'subrun', 'event', 'slc_idx']).sort_index().reset_index()
    #dfjoin = dfjoin.set_index(['run', 'subrun', 'event', 'slc_idx']).sort_index()

    #fix unit
    dfjoin["slc_opt0_time_corrected_Z_pandora"] = dfjoin["slc_opt0_time_corrected_Z_pandora"] * 1000
    
    dfsubrun = load_tree(input_file, "subruns", subrun_branches)
    dfsubrun = dfsubrun.set_index(['run','subrun']).sort_index().reset_index()
    
    dfmevprtl = load_tree(input_file, "events", mevprtl_branches)
    dfmevprtl = dfmevprtl.set_index(['run','subrun', 'event']).sort_index().reset_index()
    
    dfmct = load_tree(input_file, "events", mct_branches)
    dfmct = dfmct.set_index(['run','subrun', 'event']).sort_index().reset_index()


    if ftype == "hnl":

        scale_pot_hnl, hnl_spill = calc_scaling_pot(dfsubrun, dfjoin)

        true_counts, true_nonfv_counts = get_true_signal_in_all_spills(dfmct, scale_pot_hnl)
        start_counts, start_nonfv_counts = get_reco_signal_in_all_spills(dfjoin, scale_pot_hnl)
        
        #true reco
        print("true signals = " + str(true_counts))
        print("true nonfv signals = " + str(true_nonfv_counts))
        print("total true signals = " +str(true_counts +true_nonfv_counts))
        
        #start reco signals
        print("start signals = " + str(start_counts))
        print("start nonfv signals = " + str(start_nonfv_counts))
        print("total start signals = " +str(start_counts + start_nonfv_counts))

    elif ftype == "nu":

        dfjoin['slc_true_event_type'][dfjoin['slc_true_event_type'] == 0] = -1

    elif ftype == "cos":

        dfjoin['slc_true_event_type'][dfjoin['slc_true_event_type'] == -1] = 9

    #dfjoin = PreCut(dfjoin)
    dfjoin = cutClearCosmics(dfjoin)

    dfjoin.to_hdf(output_file, key='slc', mode = 'w')
    dfsubrun.to_hdf(output_file, key='subrun')
    dfmct.to_hdf(output_file, key='mct')
    dfmevprtl.to_hdf(output_file, key='mevprtl')

def main(args):

    start = time.time()

    load_df(args.input, args.output, args.type)
    
    now = time.time()

    print("Done converting root to dataframe, elapsed time = " + str(now -start))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()      
    parser.add_argument("--input", required=True, help="Input root file e.g. hnl.root")
    parser.add_argument("--output", required=True, help="Output hdf5 file e.g. hnl.h5")
    parser.add_argument("--type", required=True, help="hnl or nu or cos?")

    args = parser.parse_args()                             
    main(args)                                             
