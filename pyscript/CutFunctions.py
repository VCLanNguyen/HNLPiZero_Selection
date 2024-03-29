import sys
import os

os.nice(20)

import numpy as np

# Local helper script                             
hnlDIR = os.environ['_']                          
sys.path.append(hnlDIR + '/pyscript/')            
                                                
from HelperFunctions import *

#------------------------------------------------------------------------------------------------------------------#
def PreCut(df):
    #Remove clear cosmics
    df = cutClearCosmics(df)

    #Valid Opt0
    df = cutOpt0(df) 

    #Only Keep Higher Energy Shower PFP per slice
    #df = df[df['slc_pfp_shower_energy'] == df.groupby(['run','subrun','event','slc_idx'])["slc_pfp_shower_energy"].transform(max)]

    ##Valid Shower Object: Length
    #df = df[df['slc_pfp_shower_length'] > 0]
 
    ##Valid Shower Object: Energy
    #df = df[df['slc_pfp_shower_energy'] > 0]

    ##Valid Shower Object: Direction
    #df = df[(df['slc_pfp_shower_theta'] >= 0) & (df['slc_pfp_shower_theta'] <= 180)]

    ##Valid Track Object: Length
    #df = df[df['slc_pfp_track_length'] > 0]

    ##Valid Track Object: Energy
    #df = df[(df['slc_pfp_track_ke'] > 0) & (df['slc_pfp_track_ke'] < 1e9)]
    #
    ##Valid Track Object: Direction
    #df = df[(df['slc_pfp_track_theta'] >= 0) & (df['slc_pfp_track_theta'] <= 180)]

    #df = df.drop(['slc_is_clear_cosmics'], axis = 1)

    return df
#------------------------------------------------------------------------------------------------------------------#
def cutClearCosmics(df
                    , ifNotCosmics = True
                    ): 

    #1. Slice is not clear cosmics
    
    if ifNotCosmics:
        when_NotCosmics = df['slc_is_clear_cosmics'] == 0
        df = df[when_NotCosmics]

    df = df.drop(columns=['slc_is_clear_cosmics'])
    
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
def cutTrackScoreMore(df
                  , shwScore = 0.6
                  , ifShower = True
                  ):

    #3. Has at least 1 pfp with trackScore < threshold
    if ifShower:

        conditions = [
                (df['slc_pfp_track_score'] >= shwScore)
                , (df['slc_pfp_track_score'] < shwScore)
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
def cutTrackScoreLess(df
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

    df["trk_exit"] = np.where(((df["slc_pfp_track_score"] >= 0.51) & 
                                ((abs(df[trk_col+"x"]) > 195) |
                                    (df[trk_col+"y"] < -195)  | (df[trk_col+"y"] > 195) &
                                    (df[trk_col+"z"] < 5)    | (df[trk_col+"z"] > 495))),
                                    1,df["trk_exit"])

    df["shw_exit"] = np.where(((df["slc_pfp_track_score"] < 0.51) & 
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
            , muonScore = 0.035
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
            , nProton = 0
            , protonScore = 0.06
            , ifnProton = True
            , ifScore = True
            ): 
    #1. nRazzled Proton w/ threshold cut == 0
    #2. Contains no pfp with Razzled Proton Score > val
    
    if ifnProton:
        when_nProton = df['slc_n_razzled_protons_thresh'] <= nProton
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
            , nPion = 2
            , pionScore = 0.14
            , ifnPion = True
            , ifScore = True
            ): 
    #1. nRazzled Pion w/ threshold cut < 2
    #2. Contains no pfp with Razzled Pion Score > val
    
    if ifnPion:
        when_nPion = df['slc_n_razzled_pions_thresh'] <= nPion 
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
def cutElectronLess(df
            , nElectron = 2
            , electronScore = 0.14
            , ifnElectron = True
            , ifScore = True
            ): 
    #1. nRazzled Electron w/ threshold cut < 2
    #2. Contains no pfp with Razzled Electron Score > val
    
    if nElectron:
        when_nElectron = df['slc_n_razzled_electrons'] <= nElectron
        df = df[when_nElectron]
    
    if ifScore:

        conditions = [
                (df['slc_pfp_razzled_electron_score'] >= electronScore)
                , (df['slc_pfp_razzled_electron_score'] < electronScore)
                ]

        values = [1, 0]

        # create a new column and use np.select to assign values to it using our lists as arguments
        df['hasElectron'] = np.select(conditions, values)

        # creates a new dataframe that contains the column "tpc_sum"
        dfElectron = df.groupby(["run","subrun","event","slc_idx"]).agg(electron_sum = ('hasElectron','sum')).reset_index()

        # select slices that have only 
        df_noElectron = dfElectron.query("electron_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noElectron_idx = df_noElectron[["run","subrun","event","slc_idx"]]

        # merge this index with your original dataframe
        df = df.merge(df_noElectron_idx, on = ['run','subrun', 'event','slc_idx'])

        df = df.drop(columns=['hasElectron'])
    return df
#------------------------------------------------------------------------------------------------------------------#
def cutElectronMore(df
            , nElectron = 2
            , electronScore = 0.14
            , ifnElectron = True
            , ifScore = True
            ): 
    #1. nRazzled Electron w/ threshold cut < 2
    #2. Contains no pfp with Razzled Electron Score > val
    
    if nElectron:
        when_nElectron = df['slc_n_razzled_electrons'] <= nElectron
        df = df[when_nElectron]
    
    if ifScore:

        conditions = [
                (df['slc_pfp_razzled_electron_score'] <= electronScore)
                , (df['slc_pfp_razzled_electron_score'] > electronScore)
                ]

        values = [1, 0]

        # create a new column and use np.select to assign values to it using our lists as arguments
        df['hasElectron'] = np.select(conditions, values)

        # creates a new dataframe that contains the column "tpc_sum"
        dfElectron = df.groupby(["run","subrun","event","slc_idx"]).agg(electron_sum = ('hasElectron','sum')).reset_index()

        # select slices that have only 
        df_noElectron = dfElectron.query("electron_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noElectron_idx = df_noElectron[["run","subrun","event","slc_idx"]]

        # merge this index with your original dataframe
        df = df.merge(df_noElectron_idx, on = ['run','subrun', 'event','slc_idx'])

        df = df.drop(columns=['hasElectron'])
    return df
#------------------------------------------------------------------------------------------------------------------#
def cutPhotonLess(df
            , nPhoton = 2
            , photonScore = 0.14
            , ifnPhoton = True
            , ifScore = True
            ): 
    #1. nRazzled Photon w/ threshold cut < 2
    #2. Contains no pfp with Razzled Photon Score > val
    
    if nPhoton:
        when_nPhoton = df['slc_n_razzled_photons'] <= nPhoton
        df = df[when_nPhoton]
    
    if ifScore:

        conditions = [
                (df['slc_pfp_razzled_photon_score'] >= photonScore)
                , (df['slc_pfp_razzled_photon_score'] < photonScore)
                ]

        values = [1, 0]

        # create a new column and use np.select to assign values to it using our lists as arguments
        df['hasPhoton'] = np.select(conditions, values)

        # creates a new dataframe that contains the column "tpc_sum"
        dfPhoton = df.groupby(["run","subrun","event","slc_idx"]).agg(photon_sum = ('hasPhoton','sum')).reset_index()

        # select slices that have only 
        df_noPhoton = dfPhoton.query("photon_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noPhoton_idx = df_noPhoton[["run","subrun","event","slc_idx"]]

        # merge this index with your original dataframe
        df = df.merge(df_noPhoton_idx, on = ['run','subrun', 'event','slc_idx'])

        df = df.drop(columns=['hasPhoton'])
    return df
#------------------------------------------------------------------------------------------------------------------#
def cutPhotonMore(df
            , nPhoton = 2
            , photonScore = 0.14
            , ifnPhoton = True
            , ifScore = True
            ): 
    #1. nRazzled Photon w/ threshold cut < 2
    #2. Contains no pfp with Razzled Photon Score > val
    
    if nPhoton:
        when_nPhoton = df['slc_n_razzled_photons'] <= nPhoton
        df = df[when_nPhoton]
    
    if ifScore:

        conditions = [
                (df['slc_pfp_razzled_photon_score'] <= photonScore)
                , (df['slc_pfp_razzled_photon_score'] > photonScore)
                ]

        values = [1, 0]

        # create a new column and use np.select to assign values to it using our lists as arguments
        df['hasPhoton'] = np.select(conditions, values)

        # creates a new dataframe that contains the column "tpc_sum"
        dfPhoton = df.groupby(["run","subrun","event","slc_idx"]).agg(photon_sum = ('hasPhoton','sum')).reset_index()

        # select slices that have only 
        df_noPhoton = dfPhoton.query("photon_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noPhoton_idx = df_noPhoton[["run","subrun","event","slc_idx"]]

        # merge this index with your original dataframe
        df = df.merge(df_noPhoton_idx, on = ['run','subrun', 'event','slc_idx'])

        df = df.drop(columns=['hasPhoton'])
    return df
#------------------------------------------------------------------------------------------------------------------#
def cutThetaAngle(df
            , thetaAngle = 30
            , ifDir = True
            ): 
            
    if ifDir:
        conditions = [
                (df['slc_pfp_shower_theta'] <= thetaAngle)
                , (df['slc_pfp_shower_theta'] > thetaAngle)
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
def cutBeamGate(df
                , ifOpt0BeamGate = True
               ): 
        
    if ifOpt0BeamGate:
        when_Opt0Score = df['slc_opt0_time_corrected_Z_pandora'] > 367
        df = df[when_Opt0Score]

    return df
#------------------------------------------------------------------------------------------------------------------#
def cutOpt0Frac(df
                , Opt0FracMore = 0.3
                , ifOpt0FracMore = True
                , Opt0FracLess = 0.3
                , ifOpt0FracLess = True
               ): 

    if ifOpt0FracMore:
        when_Opt0Frac = df['slc_opt0_frac'] > Opt0FracMore
        df = df[when_Opt0Frac]
        
    if ifOpt0FracLess:
        when_Opt0Frac = df['slc_opt0_frac'] < Opt0FracLess
        df = df[when_Opt0Frac]

    return df
#------------------------------------------------------------------------------------------------------------------#
def cutOpt0Score(df
                , Opt0Score = 500
                , ifOpt0Valid = True
                , ifOpt0Score = True
               ): 
    
    if ifOpt0Valid:
        when_Opt0Valid = df['slc_opt0_score'] > 0 
        df = df[when_Opt0Valid]
    
    if ifOpt0Score:
        when_Opt0Valid = df['slc_opt0_score'] > Opt0Score 
        df = df[when_Opt0Valid]

    return df