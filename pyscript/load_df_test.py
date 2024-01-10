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

def load_tree(tname, branches):
   
    start = time.time() 
    folder = "hnlpizeroana"
    tree = uproot.dask("./root_files/hnl_5k.root:hnlpizeroana/events"
                       , library = "np"
                       , step_size = "100 MB"
                        , branches = branches
                       )

    now = time.time()
    print(now - start)
    test = tree['slc_crumbs_score']
    test =test.compute()
    print(test)
    now1 = time.time()
    print(now1 - now)
    #df = tree.arrays(branches, library = "pd")#, entry_start = 0, entry_stop = 5)

    return df
    
def main(args):

    branches =[
                'run'
                , 'subrun'
                , 'event'
                #, "n_slc"
                #, 'slc_primary_pfp_id'
                #, 'slc_n_pfps'
                , 'slc_crumbs_score'
                , 'slc_is_clear_cosmics'
                , 'slc_comp'
                #, 'slc_n_trks'
                #, 'slc_n_shws'
                ]

    
    dfslc = load_tree("events", branches)

    #branches =[
    #            'run'
    #            , 'subrun'
    #            , 'event'
    #            , 'slc_pfp_id'
    #            #, 'slc_pfp_track_score'
    #            , 'slc_pfp_track_start_x'
    #            , 'slc_pfp_track_start_y'
    #            , 'slc_pfp_track_start_z'
    #            #, 'slc_pfp_track_end_x'
    #            #, 'slc_pfp_track_end_y'
    #            #, 'slc_pfp_track_end_z'
    #            #, 'slc_pfp_track_dazzle_muon_score'
    #            #, 'slc_pfp_shower_razzle_photon_score'
    #            ]
    #
    #dfpfp = load_tree(args.input, "events", branches)
    #
    ##add slc idx for joining
    #dfslc["slc_idx"] = dfslc.groupby(['run', 'subrun', 'event']).transform("cumcount").add(1) - 1
    ##every row = slice
    #dfslc = dfslc.set_index(['run', 'subrun', 'event']).reset_index()
    #
    ##explode every row = slice
    #dfpfp = dfpfp.set_index(['run', 'subrun', 'event']).apply(pd.Series.explode).reset_index()
    ##change stlvector to array
    #dfpfp = dfpfp.set_index(['run', 'subrun', 'event']).applymap(lambda x: np.array(x))
    #dfpfp["slc_idx"] = dfpfp.groupby(['run', 'subrun', 'event']).transform("cumcount").add(1) - 1
    #dfpfp = dfpfp.reset_index()
   
    ##merge on slice level --> this is on purpose since we count slice, just to make life easier
    #dfjoin = pd.concat([dfslc, dfpfp], axis = 1)
    #dfjoin = dfjoin.loc[:, ~dfjoin.columns.duplicated()]
    #dfjoin = dfjoin.apply(pd.Series.explode) #explode to pfp level
    #dfjoin = dfjoin[dfjoin['slc_is_clear_cosmics'] == 0 ]
    #
    ##copy and keep a few events for testing
    #df = dfjoin
    #df = df[df['run'] == 1 ]
    #df = df[df['subrun'] == 1 ]
    #df = df[df['event'] < 4 ]
    #df = df

    ##reset index
    ##df = df.set_index(['run', 'subrun', 'event', 'slc_idx']).sort_index().reset_index()
    ##df = df.set_index(['run', 'subrun', 'event', 'slc_idx']).sort_index()
   
    ##get slice
    ##nodup_df = df.drop_duplicates(subset=["run","subrun","event", "slc_idx"])
    ##keep only relevant columns
    ##nodup_df = nodup_df[['run', 'subrun', 'event', 'slc_idx', 'slc_comp']]

    ##apply some conditions
    #conditions = [
    #            (df['slc_pfp_track_start_x'] < 0) &  (df['slc_pfp_track_start_y'] != 0) & (df['slc_pfp_track_start_z'] != 0)
    #            , (df['slc_pfp_track_start_x'] > 0) &  (df['slc_pfp_track_start_y'] != 0) & (df['slc_pfp_track_start_z'] != 0)
    #            ]

    ## create a list of the values we want to assign for each condition
    #values = [0, 1]

    ## create a new column and use np.select to assign values to it using our lists as arguments
    #df['whichTPC'] = np.select(conditions, values)

    ## creates a new dataframe that contains the column "tpc_sum"
    #pfp_tpc_sum = df.groupby(["run","subrun","event","slc_idx"]).agg(tpc_sum = ('whichTPC','sum')).reset_index()

    ## select slices that have only 
    #pfp_tpc0    = pfp_tpc_sum.query("tpc_sum==0") 
    ## alternatively you could also do (it's just longer):
    ## pfp_tpc 0 = pfp_tpc_sum[pfp_tpc_sum.tpc_sum==0]

    ## keep only the columns relevant for your index
    #pfp_tpc0_idx = pfp_tpc0[["run","subrun","event","slc_idx"]]

    ## merge this index with your original dataframe
    #tpc0_df = df.merge(pfp_tpc0_idx, on = ['run','subrun', 'event','slc_idx'])

    #print('all slice')
    #print(df)
    #print(pfp_tpc_sum)
    #print(pfp_tpc0)
    #print(pfp_tpc0_idx)
    #print(tpc0_df)
#    print(nodup_df)
#    print(nodup_df[nodup_df['slc_comp'] < 0.5])
#    print(len(nodup_df['slc_comp']))
#    print(len(nodup_df['slc_comp'][nodup_df['slc_comp'] < 0.5]))

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	args = parser.parse_args()
	main(args)
