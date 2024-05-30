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
        
    #SLICE AND PFP NEED SPECIAL HANDLING
    #print(slc_branches)
    dfslc = load_tree(input_file, "events", slc_branches)
    dfpfp = load_tree(input_file, "events", pfp_branches)
   
    #add slc idx for joining
    #dfslc["slc_idx"] = dfslc.groupby(['run', 'subrun', 'event']).transform("cumcount").add(1) - 1
    #every row = slice
    dfslc = dfslc.set_index(['run', 'subrun', 'event','slc_id']).reset_index()
   
    #explode every row = slice
    dfpfp = dfpfp.set_index(['run', 'subrun', 'event']).apply(pd.Series.explode).reset_index()
    #change stlvector to array
    dfpfp = dfpfp.set_index(['run', 'subrun', 'event']).applymap(lambda x: np.array(x))
    #dfpfp = dfpfp[~dfpfp.duplicated(['run','subrun','event','slc_id'])]
    #add slc idx for joining
    dfpfp["slc_id"] = dfpfp.groupby(['run', 'subrun', 'event']).transform("cumcount").add(1) - 1
    dfpfp = dfpfp.reset_index()
   
    #merge pfp and slc dataframe on slice idx
    dfjoin = pd.concat([dfslc, dfpfp], axis = 1)
    #dfjoin = dfslc.merge(dfpfp, how='inner', on=['run', 'subrun', 'event', 'slc_id'])

    del dfpfp

    dfjoin = dfjoin.loc[:, ~dfjoin.columns.duplicated()]
    #explode to pfp level
    dfjoin = dfjoin.apply(pd.Series.explode) #explode to pfp level
	
    dfsubrun = load_tree(input_file, "subruns", subrun_branches)
    dfmct = load_tree(input_file, "events", mct_branches)
    #dfmevprtl = load_tree(input_file, "mevprtl", mevprtl_branches)

    dfslc = dfslc[['run', 'subrun', 'event', 'slc_id', 'slc_comp']]
    dfweight = load_tree(input_file, "events", weight_branches)

    #explode every row = slice
    dfweight = dfweight.set_index(['run', 'subrun', 'event']).apply(pd.Series.explode).reset_index()
    #change stlvector to array
    dfweight = dfweight.set_index(['run', 'subrun', 'event']).applymap(lambda x: np.array(x))
    #add slc idx for joining
    dfweight["slc_id"] = dfweight.groupby(['run', 'subrun', 'event']).transform("cumcount").add(1) - 1
    dfweight = dfweight.reset_index()

    #merge weight and slc dataframe on slice idx
    dfweight = pd.concat([dfslc, dfweight], axis = 1)
    dfweight = dfweight.loc[:, ~dfweight.columns.duplicated()]

    del dfslc

    return dfjoin, dfsubrun, dfmct, dfweight

def manipulate_dfjoin(record, dfjoin, dfsubrun, dfmct, ftype):

    #fix unit
    dfjoin["slc_opt0_time_corrected_Z_pandora"] = dfjoin["slc_opt0_time_corrected_Z_pandora"] * 1000

    if ftype == "hnl":

        #scale_pot_hnl, hnl_spill = calc_scaling_pot(dfsubrun, dfjoin)
        true_counts, true_nonfv_counts = get_true_signal_in_all_spills(dfmct, 1)
        start_counts, start_nonfv_counts = get_reco_signal_in_all_spills(dfjoin, 1)
        
        #true reco
        record.write("true_signals = " + str(true_counts) + "\n")
        record.write("true_nonfv_signals = " + str(true_nonfv_counts) + "\n")
        record.write("total_true_signals = " +str(true_counts +true_nonfv_counts) + "\n")
        
        #start reco signals
        record.write("start_signals = " + str(start_counts) + "\n")
        record.write("start_nonfv_signals = " + str(start_nonfv_counts) + "\n")
        record.write("total_start_signals = " +str(start_counts + start_nonfv_counts) + "\n")

    elif ftype == "nu":
        dfjoin['slc_true_event_type'][dfjoin['slc_true_event_type'] == 0] = -1

    elif ftype == "cos":
        dfjoin['slc_true_event_type'][dfjoin['slc_true_event_type'] == -1] = 9

    #Apply cutting on slice dataframe
    dfjoin = cutClearCosmics(dfjoin)
    
    return dfjoin

def save_hdf5(dfjoin, dfsubrun, dfmct, output_file):

    dfjoin.to_hdf(output_file, key='slc', mode = 'w')
    dfsubrun.to_hdf(output_file, key='subrun')
    dfmct.to_hdf(output_file, key='mct')
    #dfmevprtl.to_hdf(output_file, key='mevprtl')

def make_new_dataframe():
    
    dfjoin_cc  = pd.DataFrame() 
    dfsubrun_cc  = pd.DataFrame() 
    dfmct_cc  = pd.DataFrame() 
    #dfmevprtl_cc  = pd.DataFrame()
    dfweight_cc = pd.DataFrame()

    return dfjoin_cc, dfsubrun_cc, dfmct_cc, dfweight_cc

def reset_index_dataframe(dfjoin_cc, dfsubrun_cc, dfmct_cc, dfweight_cc):
    ##reset index
    dfjoin_cc = dfjoin_cc.set_index(['run','subrun', 'event']).sort_index().reset_index()
    dfsubrun_cc = dfsubrun_cc.set_index(['run','subrun']).sort_index().reset_index()
    dfmct_cc = dfmct_cc.set_index(['run','subrun', 'event']).sort_index().reset_index()
    dfweight_cc = dfweight_cc.set_index(['run','subrun', 'event']).sort_index().reset_index()

    return dfjoin_cc, dfsubrun_cc, dfmct_cc, dfweight_cc

def main(args):

    file_list = open(args.input, "r")
    lines = file_list.readlines()
    file_list.close()
    nlines = len(lines) - 1
    
    start = time.time()
    
    dfjoin_cc, dfsubrun_cc, dfmct_cc, dfweight_cc = make_new_dataframe()
    itr = 0
    
    record = open(args.output+"_slc_count.txt", "w")
    record.write("Save Slice Counts Before Clear Cosmics Cut \n")

    for i, line in enumerate(lines):

        print("Reading in files " + str(i) + ": " + str(line))
        dfjoin, dfsubrun, dfmct, dfweight = load_df(line)
   
        dfjoin_cc = pd.concat([dfjoin_cc, dfjoin], ignore_index = True) 
        dfsubrun_cc = pd.concat([dfsubrun_cc, dfsubrun], ignore_index = True) 
        dfmct_cc = pd.concat([dfmct_cc, dfmct], ignore_index = True) 
        dfweight_cc = pd.concat([dfweight_cc, dfweight], ignore_index = True) 

        if ((i != 0) and (i%20 == 0)) or ((i != 0) and (i%nlines == 0)) :
    
            print ("Saving iteration " + str(itr) + ", file counter " + str(i))
            
            reset_index_dataframe(dfjoin_cc, dfsubrun_cc, dfmct_cc, dfweight_cc)       
    
            ##fix stuff in here
            record.write("Iteration " + str(itr) + " \n")
            dfjoin_cc = manipulate_dfjoin(record, dfjoin_cc, dfsubrun_cc, dfmct_cc, args.type)

            ##save to pickle
            dfjoin_cc.to_pickle(args.output+"_slc_" + str(itr) + ".pkl", protocol = 5)
            dfsubrun_cc.to_pickle(args.output+"_subrun_" + str(itr) + ".pkl", protocol = 5)
            dfmct_cc.to_pickle(args.output+"_mct_" + str(itr) + ".pkl", protocol = 5)
            dfweight_cc.to_pickle(args.output+"_flxw_" + str(itr) + ".pkl", protocol = 5)
    
            del dfjoin_cc
            del dfsubrun_cc
            del dfmct_cc
            del dfweight_cc
            
            ##make new dataframe and restart the cycle
            dfjoin_cc, dfsubrun_cc, dfmct_cc, dfweight_cc = make_new_dataframe()
            itr += 1

    record.close()

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
