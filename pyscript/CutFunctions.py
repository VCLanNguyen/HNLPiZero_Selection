import sys
import os

os.nice(20)

import numpy as np

# Local helper script                             
hnlDIR = os.environ['_']                          
sys.path.append(hnlDIR + '/pyscript/')            
                                                
from HelperFunctions import *
#------------------------------------------------------------------------------------------------------------------#
def cutClearCosmics(df
                    , ifNotCosmics = True
                    ): 

    #1. Slice is not clear cosmics
    
    if ifNotCosmics:
        when_NotCosmics = df['slc_is_clear_cosmics'] == 0
        df = df[when_NotCosmics]

    return df
    
#------------------------------------------------------------------------------------------------------------------#

def cutFV(df
          , ifRecoVtxFV = True
          ): 
    
    #2. Reco vertex is contained in FV
    if ifRecoVtxFV:
        when_RecoVtxFV = df['slc_is_fv'] == 1
        df = df[when_RecoVtxFV]

    return df
    
#------------------------------------------------------------------------------------------------------------------#
def cutTrackScore(df
                  , shwScore = 0.6
                  , ifShower = True
                  ):

    #3. Has at least 1 pfp with trackScore < threshold
    if ifShower:

        conditions = [
                (df['slc_pfp_track_score'] <= shwScore)
                , (df['slc_pfp_track_score'] > shwScore)
                ]

        values = [1, 0]

        # create a new column and use np.select to assign values to it using our lists as arguments
        df['showerLike'] = np.select(conditions, values)

        # creates a new dataframe that contains the column "tpc_sum"
        dfShower = df.groupby(["run","subrun","event","slc_idx"]).agg(shower_sum = ('showerLike','sum')).reset_index()

        # select slices that have only 
        dfShower = dfShower.query("shower_sum > 0") 
        
        # keep only the columns relevant for your index
        dfShower_idx = dfShower[["run","subrun","event","slc_idx"]]

        # merge this index with your original dataframe
        df = df.merge(dfShower_idx, on = ['run','subrun', 'event','slc_idx'])
        
        #drop useless columns
        df = df.drop(columns=['showerLike'])

    return df
    
#------------------------------------------------------------------------------------------------------------------#
def cutCosmics(df
                , crumbsScore = 0.00
                , ifCrumbs = True
               ): 

    #1. CRUMBS Score > threshold
    
    if ifCrumbs:
        when_Crumbs = df['slc_crumbs_score'] > crumbsScore 
        df = df[when_Crumbs]

    return df

#------------------------------------------------------------------------------------------------------------------#
def cutContainment(df
                   , ifTrkCont = True
                   , ifShwCont = True
                   , trk_col ="slc_pfp_track_end_"
                   , shw_col ="slc_pfp_shower_end_"
                   ):

    #1. Containment: no end point within 5 cm of detector

    df["shw_exit"] = 0 
    df["trk_exit"] = 0

    df["trk_exit"] = np.where(((df["slc_pfp_track_score"] >= 0.5) & 
                                ((abs(df[trk_col+"x"]) > 195) |
                                    (df[trk_col+"y"] < -195)  | (df[trk_col+"y"] > 195) &
                                    (df[trk_col+"z"] < 5)    | (df[trk_col+"z"] > 495))),
                                    1,df["trk_exit"])

    df["shw_exit"] = np.where(((df["slc_pfp_track_score"] < 0.5) & 
                            ((abs(df[shw_col+"x"]) > 195) |
                                (df[shw_col+"y"] < -195)  | (df[shw_col+"y"] > 195) &
                                (df[shw_col+"z"] < 5)    | (df[shw_col+"z"] > 495))),
                                1,df["shw_exit"])

    # sum the number of exiting trks/shws 
    df_exit = df.groupby(["run","subrun","event", "slc_idx"]).agg(ntrk_exit = ('trk_exit','sum'),
                                                                      nshw_exit = ('shw_exit','sum')).reset_index()

    # require that there are no exiting trks/shw if specified
    if ifTrkCont: 
        df_trk_cont = df_exit.query("ntrk_exit==0")
        
        # keep only the columns relevant for your index
        df_trk_cont_idx = df_trk_cont[["run","subrun","event","slc_idx"]]

        # merge this index with your original dataframe
        df = df.merge(df_trk_cont_idx, on = ['run','subrun', 'event','slc_idx'])

    # require that there are no exiting trks/shw if specified
    if ifShwCont: 
        df_shw_cont = df_exit.query("nshw_exit==0")
        
        # keep only the columns relevant for your index
        df_shw_cont_idx = df_shw_cont[["run","subrun","event","slc_idx"]]

        # merge this index with your original dataframe
        df = df.merge(df_shw_cont_idx, on = ['run','subrun', 'event','slc_idx'])
    
    #drop useless columns
    df = df.drop(columns=['shw_exit', 'trk_exit'])
   
    return df

#------------------------------------------------------------------------------------------------------------------#
def cutMuon(df
            , muonScore = 0.04
            , ifnMuon = True
            , ifScore = True
            ): 

    #1. nRazzled Muon = 0
    #2. contains no pfp with Razzled Muon Score > val

    if ifnMuon:
        when_nMuon = df['slc_n_razzled_muons'] == 0 
        df = df[when_nMuon]
    
    if ifScore:

        conditions = [
                (df['slc_pfp_razzled_muon_score'] >= muonScore)
                , (df['slc_pfp_razzled_muon_score'] < muonScore)
                ]

        values = [1, 0]

        # create a new column and use np.select to assign values to it using our lists as arguments
        df['hasMuon'] = np.select(conditions, values)

        # creates a new dataframe that contains the column "tpc_sum"
        dfMuon = df.groupby(["run","subrun","event","slc_idx"]).agg(muon_sum = ('hasMuon','sum')).reset_index()

        # select slices that have only 
        df_noMuon = dfMuon.query("muon_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noMuon_idx = df_noMuon[["run","subrun","event","slc_idx"]]

        # merge this index with your original dataframe
        df = df.merge(df_noMuon_idx, on = ['run','subrun', 'event','slc_idx'])
        
        df = df.drop(columns=['hasMuon'])

    return df

#------------------------------------------------------------------------------------------------------------------#
def cutProton(df
            , protonScore = 0.06
            , ifnProton = True
            , ifScore = True
            ): 
    #1. nRazzled Proton w/ threshold cut == 0
    #2. Contains no pfp with Razzled Proton Score > val
    
    if ifnProton:
        when_nProton = df['slc_n_razzled_protons_thresh'] == 0 
        df = df[when_nProton]
    
    if ifScore:

        conditions = [
                (df['slc_pfp_razzled_proton_score'] >= protonScore)
                , (df['slc_pfp_razzled_proton_score'] < protonScore)
                ]

        values = [1, 0]

        # create a new column and use np.select to assign values to it using our lists as arguments
        df['hasProton'] = np.select(conditions, values)

        # creates a new dataframe that contains the column "tpc_sum"
        dfProton = df.groupby(["run","subrun","event","slc_idx"]).agg(proton_sum = ('hasProton','sum')).reset_index()

        # select slices that have only 
        df_noProton = dfProton.query("proton_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noProton_idx = df_noProton[["run","subrun","event","slc_idx"]]

        # merge this index with your original dataframe
        df = df.merge(df_noProton_idx, on = ['run','subrun', 'event','slc_idx'])
        
        df = df.drop(columns=['hasProton'])

    return df
#------------------------------------------------------------------------------------------------------------------#
def cutPion(df
            , pionScore = 0.14
            , ifnPion = True
            , ifScore = True
            ): 
    #1. nRazzled Pion w/ threshold cut < 2
    #2. Contains no pfp with Razzled Pion Score > val
    
    if ifnPion:
        when_nPion = df['slc_n_razzled_pions_thresh'] == 0 
        df = df[when_nPion]
    
    if ifScore:

        conditions = [
                (df['slc_pfp_razzled_pion_score'] >= pionScore)
                , (df['slc_pfp_razzled_pion_score'] < pionScore)
                ]

        values = [1, 0]

        # create a new column and use np.select to assign values to it using our lists as arguments
        df['hasPion'] = np.select(conditions, values)

        # creates a new dataframe that contains the column "tpc_sum"
        dfPion = df.groupby(["run","subrun","event","slc_idx"]).agg(pion_sum = ('hasPion','sum')).reset_index()

        # select slices that have only 
        df_noPion = dfPion.query("pion_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noPion_idx = df_noPion[["run","subrun","event","slc_idx"]]

        # merge this index with your original dataframe
        df = df.merge(df_noPion_idx, on = ['run','subrun', 'event','slc_idx'])

        df = df.drop(columns=['hasPion'])
    return df
#------------------------------------------------------------------------------------------------------------------#
def cutShower(df
            , angle = 30
            , ifnShower = True
            , ifDir = True
            ): 

    #1. Contains at least 1 Razzled photon
    #2. Contains at least 1 shower with opening angle close 2 beams
    print(df['slc_n_razzled_photons'])

    if ifnShower:
        when_nShower = df['slc_n_razzled_photons'] > 0 
        df = df[when_nShower]
    print(df['slc_n_razzled_photons'])
            
    if ifDir:
        conditions = [
                (df['slc_pfp_angle2Beam'] <= angle)
                , (df['slc_pfp_angle2Beam'] > angle)
                ]

        values = [1, 0]

        # create a new column and use np.select to assign values to it using our lists as arguments
        df['close2Beam'] = np.select(conditions, values)

        # creates a new dataframe that contains the column "tpc_sum"
        dfBeam = df.groupby(["run","subrun","event","slc_idx"]).agg(beam_sum = ('close2Beam','sum')).reset_index()

        # select slices that have only 
        dfBeam = dfBeam.query("beam_sum > 0") 
        
        # keep only the columns relevant for your index
        dfBeam_idx = dfBeam[["run","subrun","event","slc_idx"]]

        # merge this index with your original dataframe
        df = df.merge(dfBeam_idx, on = ['run','subrun', 'event','slc_idx'])

        df = df.drop(columns=['close2Beam'])

    return df
#------------------------------------------------------------------------------------------------------------------#

def cutBetweenBucket(df
                     , ifCut = True
                    ):
    if ifCut:
        lb_arr, ub_arr = make_interval(375, 380)
        #lb_arr, ub_arr = make_interval(377, 380)
        
        df['isBetweenBucket'] = df.apply(lambda row : checkInterval(row['slc_opt0_time_corrected_Z_pandora'], lb_arr, ub_arr), axis = 1)
        when_isBetweenBucket = df['isBetweenBucket'] == True

        df = df[when_isBetweenBucket]
        
        df = df.drop(columns=['isBetweenBucket'])
   
    return df
#------------------------------------------------------------------------------------------------------------------#
def cutOpt0Score(df
                , Opt0Score = 500
                , ifOpt0Score = True
               ): 

    #1. Opt0 Score > threshold
    
    if ifOpt0Score:
        when_Opt0Score = df['slc_opt0_score'] > Opt0Score 
        df = df[when_Opt0Score]

    return df
