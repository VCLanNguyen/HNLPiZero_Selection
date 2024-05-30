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
    #df = df[df['slc_pfp_shower_energy'] == df.groupby(['run','subrun','event','slc_id'])["slc_pfp_shower_energy"].transform(max)]

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
        dfShower = df.groupby(["run","subrun","event","slc_id"]).agg(shower_sum = ('showerLike','sum')).reset_index()

        # select slices that have only 
        dfShower = dfShower.query("shower_sum > 0") 
        
        # keep only the columns relevant for your index
        dfShower_idx = dfShower[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(dfShower_idx, on = ['run','subrun', 'event','slc_id'])
        
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
        dfShower = df.groupby(["run","subrun","event","slc_id"]).agg(shower_sum = ('showerLike','sum')).reset_index()

        # select slices that have only 
        dfShower = dfShower.query("shower_sum > 0") 
        
        # keep only the columns relevant for your index
        dfShower_idx = dfShower[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(dfShower_idx, on = ['run','subrun', 'event','slc_id'])
        
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
def cutContainment_v2(df
                   , ifTrkCont = True
                   , ifShwCont = True
                   ):

    #1. Containment: no end point within 5 cm of detector

    df["shw_exit"] = 0 
    df["trk_exit"] = 0

    df["trk_exit"] = np.where(((df["slc_pfp_track_score"] >= 0.5) & (df["slc_pfp_track_contained"] == False)) ,
                                    1,df["trk_exit"])

    df["shw_exit"] = np.where(((df["slc_pfp_track_score"] < 0.5) & (df["slc_pfp_shower_contained"] == False)),
                                1,df["shw_exit"])

    # sum the number of exiting trks/shws 
    df_exit = df.groupby(["run","subrun","event", "slc_id"]).agg(ntrk_exit = ('trk_exit','sum'),
                                                                      nshw_exit = ('shw_exit','sum')).reset_index()

    # require that there are no exiting trks/shw if specified
    if ifTrkCont: 
        df_trk_cont = df_exit.query("ntrk_exit==0")
        
        # keep only the columns relevant for your index
        df_trk_cont_idx = df_trk_cont[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(df_trk_cont_idx, on = ['run','subrun', 'event','slc_id'])

    # require that there are no exiting trks/shw if specified
    if ifShwCont: 
        df_shw_cont = df_exit.query("nshw_exit==0")
        
        # keep only the columns relevant for your index
        df_shw_cont_idx = df_shw_cont[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(df_shw_cont_idx, on = ['run','subrun', 'event','slc_id'])
    
    #drop useless columns
    df = df.drop(columns=['shw_exit', 'trk_exit'])
   
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
    df_exit = df.groupby(["run","subrun","event", "slc_id"]).agg(ntrk_exit = ('trk_exit','sum'),
                                                                      nshw_exit = ('shw_exit','sum')).reset_index()

    # require that there are no exiting trks/shw if specified
    if ifTrkCont: 
        df_trk_cont = df_exit.query("ntrk_exit==0")
        
        # keep only the columns relevant for your index
        df_trk_cont_idx = df_trk_cont[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(df_trk_cont_idx, on = ['run','subrun', 'event','slc_id'])

    # require that there are no exiting trks/shw if specified
    if ifShwCont: 
        df_shw_cont = df_exit.query("nshw_exit==0")
        
        # keep only the columns relevant for your index
        df_shw_cont_idx = df_shw_cont[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(df_shw_cont_idx, on = ['run','subrun', 'event','slc_id'])
    
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
        dfMuon = df.groupby(["run","subrun","event","slc_id"]).agg(muon_sum = ('hasMuon','sum')).reset_index()

        # select slices that have only 
        df_noMuon = dfMuon.query("muon_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noMuon_idx = df_noMuon[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(df_noMuon_idx, on = ['run','subrun', 'event','slc_id'])
        
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
        dfProton = df.groupby(["run","subrun","event","slc_id"]).agg(proton_sum = ('hasProton','sum')).reset_index()

        # select slices that have only 
        df_noProton = dfProton.query("proton_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noProton_idx = df_noProton[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(df_noProton_idx, on = ['run','subrun', 'event','slc_id'])
        
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
        dfPion = df.groupby(["run","subrun","event","slc_id"]).agg(pion_sum = ('hasPion','sum')).reset_index()

        # select slices that have only 
        df_noPion = dfPion.query("pion_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noPion_idx = df_noPion[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(df_noPion_idx, on = ['run','subrun', 'event','slc_id'])

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
        dfElectron = df.groupby(["run","subrun","event","slc_id"]).agg(electron_sum = ('hasElectron','sum')).reset_index()

        # select slices that have only 
        df_noElectron = dfElectron.query("electron_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noElectron_idx = df_noElectron[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(df_noElectron_idx, on = ['run','subrun', 'event','slc_id'])

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
        dfElectron = df.groupby(["run","subrun","event","slc_id"]).agg(electron_sum = ('hasElectron','sum')).reset_index()

        # select slices that have only 
        df_noElectron = dfElectron.query("electron_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noElectron_idx = df_noElectron[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(df_noElectron_idx, on = ['run','subrun', 'event','slc_id'])

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
        dfPhoton = df.groupby(["run","subrun","event","slc_id"]).agg(photon_sum = ('hasPhoton','sum')).reset_index()

        # select slices that have only 
        df_noPhoton = dfPhoton.query("photon_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noPhoton_idx = df_noPhoton[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(df_noPhoton_idx, on = ['run','subrun', 'event','slc_id'])

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
        dfPhoton = df.groupby(["run","subrun","event","slc_id"]).agg(photon_sum = ('hasPhoton','sum')).reset_index()

        # select slices that have only 
        df_noPhoton = dfPhoton.query("photon_sum == 0") 
        
        # keep only the columns relevant for your index
        df_noPhoton_idx = df_noPhoton[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(df_noPhoton_idx, on = ['run','subrun', 'event','slc_id'])

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
        dfBeam = df.groupby(["run","subrun","event","slc_id"]).agg(beam_sum = ('close2Beam','sum')).reset_index()

        # select slices that have only 
        dfBeam = dfBeam.query("beam_sum > 0") 
        
        # keep only the columns relevant for your index
        dfBeam_idx = dfBeam[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(dfBeam_idx, on = ['run','subrun', 'event','slc_id'])

        df = df.drop(columns=['close2Beam'])

    return df
#------------------------------------------------------------------------------------------------------------------#
def cutShowerEnergy(df
            , energy = 500
            , ifCut = True
            ): 
            
    if ifCut:
        conditions = [
                (df['slc_pfp_shower_energy'] >= energy)
                , (df['slc_pfp_shower_energy'] < energy)
                ]

        values = [1, 0]

        # create a new column and use np.select to assign values to it using our lists as arguments
        df['highEnergy'] = np.select(conditions, values)

        # creates a new dataframe that contains the column "tpc_sum"
        dfBeam = df.groupby(["run","subrun","event","slc_id"]).agg(E_sum = ('highEnergy','sum')).reset_index()

        # select slices that have only 
        dfBeam = dfBeam.query("E_sum > 0") 
        
        # keep only the columns relevant for your index
        dfBeam_idx = dfBeam[["run","subrun","event","slc_id"]]

        # merge this index with your original dataframe
        df = df.merge(dfBeam_idx, on = ['run','subrun', 'event','slc_id'])

        df = df.drop(columns=['highEnergy'])

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
                , gateBegin = 367
                , gateEnd = 1967
               ): 
        
    if ifOpt0BeamGate:
        whenBegin = df['slc_opt0_time'] > gateBegin
        whenEnd = df['slc_opt0_time'] < gateEnd
        df = df[whenBegin & whenEnd]

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
#------------------------------------------------------------------------------------------------------------------#
def cutModt(df, score):
    when_modt = df["mod_t"] >= score
    df = df[when_modt]
    return df    
#------------------------------------------------------------------------------------------------------------------#
def vary_crumbsScore(df_hnl, df_nu, df_cos, true_counts, start_counts):
    
    step = 0.02
    cutStep = np.arange(-1.0, 1 + step, step)
    
    p_arr, eff_arr, peff_arr = [], [], []
    
    purity_start, eff_start = calc_purity_eff(df_hnl, df_nu, df_cos, true_counts, start_counts)
    
    for c in cutStep:
        df_hnl_cut = cutCosmics(df_hnl, crumbsScore = c)
        df_nu_cut = cutCosmics(df_nu, crumbsScore = c)
        df_cos_cut = cutCosmics(df_cos, crumbsScore = c)
        
        purity, eff = calc_purity_eff(df_hnl_cut, df_nu_cut, df_cos_cut, true_counts, start_counts)
        
        p_arr.append(purity)
        eff_arr.append(eff)
        peff_arr.append(purity*eff)
    
    plot_purity_eff(p_arr, eff_arr, peff_arr, cutStep, purity_start, eff_start, loc = 'center left')
#------------------------------------------------------------------------------------------------------------------#
def vary_OpT0Score(df_hnl, df_nu, df_cos, true_counts, start_counts):
    
    step = 50
    cutStep = np.arange(0, 2000 + step, step)
    
    p_arr, eff_arr, peff_arr = [], [], []
    
    purity_start, eff_start = calc_purity_eff(df_hnl, df_nu, df_cos, true_counts, start_counts)
  
    for c in cutStep:
        df_hnl_cut = cutOpt0Score(df_hnl, Opt0Score = c)
        df_nu_cut = cutOpt0Score(df_nu, Opt0Score = c)
        df_cos_cut = cutOpt0Score(df_cos, Opt0Score = c)
        
        purity, eff = calc_purity_eff(df_hnl_cut, df_nu_cut, df_cos_cut, true_counts, start_counts)
        
        p_arr.append(purity)
        eff_arr.append(eff)
        peff_arr.append(purity*eff)
     
    plot_purity_eff(p_arr, eff_arr, peff_arr, cutStep, purity_start, eff_start, loc = 'center left')
#------------------------------------------------------------------------------------------------------------------#
def vary_OpT0FracMore(df_hnl, df_nu, df_cos, true_counts, start_counts):
    
    step = 0.02
    cutStep = np.arange(-1, 0.5 + step, step)
    
    p_arr, eff_arr, peff_arr = [], [], []
    
    purity_start, eff_start = calc_purity_eff(df_hnl, df_nu, df_cos, true_counts, start_counts)
    
    for c in cutStep:
        df_hnl_cut = cutOpt0Frac(df_hnl, Opt0FracMore = c, Opt0FracLess = 999)
        df_nu_cut = cutOpt0Frac(df_nu, Opt0FracMore = c, Opt0FracLess = 999)
        df_cos_cut = cutOpt0Frac(df_cos, Opt0FracMore = c, Opt0FracLess = 999)
        
        purity, eff = calc_purity_eff(df_hnl_cut, df_nu_cut, df_cos_cut, true_counts, start_counts)
    
        p_arr.append(purity)
        eff_arr.append(eff)
        peff_arr.append(purity*eff)
     
    plot_purity_eff(p_arr, eff_arr, peff_arr, cutStep, purity_start, eff_start, loc = 'center left')
#------------------------------------------------------------------------------------------------------------------#
def vary_OpT0FracLess(df_hnl, df_nu, df_cos, true_counts, start_counts):
    
    step = 0.02
    cutStep = np.arange(0, 2 + step, step)
    
    p_arr, eff_arr, peff_arr = [], [], []
    
    purity_start, eff_start = calc_purity_eff(df_hnl, df_nu, df_cos, true_counts, start_counts)
    
    for c in cutStep:
        df_hnl_cut = cutOpt0Frac(df_hnl, Opt0FracMore = -999, Opt0FracLess = c)
        df_nu_cut = cutOpt0Frac(df_nu, Opt0FracMore = -999, Opt0FracLess = c)
        df_cos_cut = cutOpt0Frac(df_cos, Opt0FracMore = -999, Opt0FracLess = c)
        
        purity, eff = calc_purity_eff(df_hnl_cut, df_nu_cut, df_cos_cut, true_counts, start_counts)
    
        p_arr.append(purity)
        eff_arr.append(eff)
        peff_arr.append(purity*eff)
     
    plot_purity_eff(p_arr, eff_arr, peff_arr, cutStep, purity_start, eff_start, loc = 'center right')
#------------------------------------------------------------------------------------------------------------------#
def vary_MuonScore(df_hnl, df_nu, df_cos, true_counts, start_counts):
    
    step = 0.02
    cutStep = np.arange(0, 1 + step, step)
    
    p_arr, eff_arr, peff_arr = [], [], []
    
    purity_start, eff_start = calc_purity_eff(df_hnl, df_nu, df_cos, true_counts, start_counts)
    
    for c in cutStep:
        df_hnl_cut = cutMuon(df_hnl, muonScore = c)
        df_nu_cut = cutMuon(df_nu, muonScore = c)
        df_cos_cut = cutMuon(df_cos, muonScore = c)
        
        purity, eff = calc_purity_eff(df_hnl_cut, df_nu_cut, df_cos_cut, true_counts, start_counts)
        
        p_arr.append(purity)
        eff_arr.append(eff)
        peff_arr.append(purity*eff)
     
    plot_purity_eff(p_arr, eff_arr, peff_arr, cutStep, purity_start, eff_start , loc = 'center right')
#------------------------------------------------------------------------------------------------------------------#
def vary_ProtonScore(df_hnl, df_nu, df_cos, true_counts, start_counts):
    
    step = 0.02
    cutStep = np.arange(0, 1 + step, step)
    
    p_arr, eff_arr, peff_arr = [], [], []
    
    purity_start, eff_start  = calc_purity_eff(df_hnl, df_nu, df_cos, true_counts, start_counts)
    
    for c in cutStep:
        df_hnl_cut = cutProton(df_hnl, nProton = 0, protonScore = c)
        df_nu_cut = cutProton(df_nu, nProton = 0, protonScore = c)
        df_cos_cut = cutProton(df_cos, nProton = 0, protonScore = c)
        
        purity, eff = calc_purity_eff(df_hnl_cut, df_nu_cut, df_cos_cut, true_counts, start_counts)
        
        p_arr.append(purity)
        eff_arr.append(eff)
        peff_arr.append(purity*eff)
     
    plot_purity_eff(p_arr, eff_arr, peff_arr, cutStep, purity_start, eff_start, loc = 'center right')
#------------------------------------------------------------------------------------------------------------------#
def vary_PionScore(df_hnl, df_nu, df_cos, true_counts, start_counts):
    
    step = 0.02
    cutStep = np.arange(0, 1 + step, step)
    
    p_arr, eff_arr, peff_arr = [], [], []
    
    purity_start, eff_start  = calc_purity_eff(df_hnl, df_nu, df_cos, true_counts, start_counts)
    
    for c in cutStep:
        df_hnl_cut = cutPion(df_hnl, nPion = 0, pionScore = c)
        df_nu_cut = cutPion(df_nu, nPion = 0, pionScore = c)
        df_cos_cut = cutPion(df_cos, nPion = 0, pionScore = c)
        
        purity, eff = calc_purity_eff(df_hnl_cut, df_nu_cut, df_cos_cut, true_counts, start_counts)
        
        p_arr.append(purity)
        eff_arr.append(eff)
        peff_arr.append(purity*eff)
     
    plot_purity_eff(p_arr, eff_arr, peff_arr, cutStep, purity_start, eff_start, loc = 'center right')
#------------------------------------------------------------------------------------------------------------------#    
def vary_Theta(df_hnl, df_nu, df_cos, true_counts, start_counts):
    
    step = 0.2
    cutStep = np.arange(6, 50 + step, step)
    
    p_arr, eff_arr, peff_arr = [], [], []
    
    purity_start, eff_start = calc_purity_eff(df_hnl, df_nu, df_cos, true_counts, start_counts)
    
    for c in cutStep:
        df_hnl_cut = cutThetaAngle(df_hnl, thetaAngle = c) 
        df_nu_cut = cutThetaAngle(df_nu, thetaAngle = c) 
        df_cos_cut = cutThetaAngle(df_cos, thetaAngle = c) 
        
        purity, eff = calc_purity_eff(df_hnl_cut, df_nu_cut, df_cos_cut, true_counts, start_counts)
        
        p_arr.append(purity)
        eff_arr.append(eff)
        peff_arr.append(purity*eff)
    
    plot_purity_eff(p_arr, eff_arr, peff_arr, cutStep, purity_start, eff_start, loc = 'center right')
#------------------------------------------------------------------------------------------------------------------#     
def vary_Modt(df_hnl, df_nu, df_cos, true_counts, start_counts):
    
    step = 1
    cutStep = np.arange(0, 19 + step, step)
    
    p_arr, eff_arr, peff_arr = [], [], []
    
    purity_start, eff_start = calc_purity_eff(df_hnl, df_nu, df_cos, true_counts, start_counts)
    
    for c in cutStep:
        df_hnl_cut = cutModt(df_hnl, score = c) 
        df_nu_cut = cutModt(df_nu, score = c) 
        df_cos_cut = cutModt(df_cos, score = c) 
        
        purity, eff = calc_purity_eff(df_hnl_cut, df_nu_cut, df_cos_cut, true_counts, start_counts)
        
        p_arr.append(purity)
        eff_arr.append(eff)
        peff_arr.append(purity*eff)
    
    plot_purity_eff(p_arr, eff_arr, peff_arr, cutStep, purity_start, eff_start, loc = 'center right')
#------------------------------------------------------------------------------------------------------------------#     
def vary_ShowerEnergy(df_hnl, df_nu, df_cos, true_counts, start_counts):
    
    step = 50
    cutStep = np.arange(0, 1000 + step, step)
    
    p_arr, eff_arr, peff_arr = [], [], []
    
    purity_start, eff_start = calc_purity_eff(df_hnl, df_nu, df_cos, true_counts, start_counts)
    
    for c in cutStep:
        df_hnl_cut = cutShowerEnergy(df_hnl, energy = c) 
        df_nu_cut = cutShowerEnergy(df_nu, energy = c) 
        df_cos_cut = cutShowerEnergy(df_cos, energy = c) 
        
        purity, eff = calc_purity_eff(df_hnl_cut, df_nu_cut, df_cos_cut, true_counts, start_counts)
        
        p_arr.append(purity)
        eff_arr.append(eff)
        peff_arr.append(purity*eff)
       
    plot_purity_eff(p_arr, eff_arr, peff_arr, cutStep, purity_start, eff_start, loc = 'center right')
#------------------------------------------------------------------------------------------------------------------#       
def pi0mass(x):
    
    mBest = -999
    diff = 99999;
    
    for i in range(0, x['n_pfp']):
        for j in range(0, x['n_pfp']):
            if i == j:
                continue
            else:
                v1 = [x['slc_pfp_shower_dir_x'][0], x['slc_pfp_shower_dir_y'][0], x['slc_pfp_shower_dir_z'][0]]
                v2 = [x['slc_pfp_shower_dir_x'][1], x['slc_pfp_shower_dir_y'][1], x['slc_pfp_shower_dir_z'][1]]
    
                v1_u = unit_vector(v1)
                v2_u = unit_vector(v2)
    
                cosTheta = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)

                mass = np.sqrt(2 * x['slc_pfp_shower_energy'][0] * x['slc_pfp_shower_energy'][1] * (1 - cosTheta ));
                
                if abs(mass - 135) < diff:
                    mBest = mass
    return mBest
#------------------------------------------------------------------------------------------------------------------#   
def get_primary_shw_df(df):
    df['highest_energy'] = df.groupby(['run','subrun','event','slc_id'])['slc_pfp_shower_energy'].transform(max)
    df_prim = df[df['slc_pfp_shower_energy'] >= df['highest_energy']]
    
    df = df.drop(['highest_energy'], axis=1)
    df_prim = df_prim.drop(['highest_energy'], axis=1)
    
    return df_prim
#------------------------------------------------------------------------------------------------------------------#   
def apply_shower_cut(df, dfshw1, dfshw2):
    temp1 = dfshw1[['run','subrun','event','slc_id']]
    temp1 = temp1.merge(df, how='inner', on=['run','subrun','event','slc_id'])
    
    temp2 = dfshw2[['run','subrun','event','slc_id']]
    temp2 = temp2.merge(df, how='inner', on=['run','subrun','event','slc_id'])
    
    concat = pd.concat([temp1, temp2])
    concat = concat.drop_duplicates()
    
    return concat
#------------------------------------------------------------------------------------------------------------------#   
def split_my_df(df):
    
    temp = df[['run', 'subrun', 'event', 'slc_id'
                   , 'slc_pfp_shower_energy', 'slc_comp', 'scale_pot', 'slc_true_event_type', 'mod_t'
                   ,'slc_opt0_frac' ,'slc_opt0_measPE'
                   , 'slc_pfp_shower_dir_x' , 'slc_pfp_shower_dir_y' , 'slc_pfp_shower_dir_z'
                   , 'slc_pfp_shower_theta', 'slc_pfp_shower_phi'
                   , 'slc_pfp_shower_conv_gap', 'slc_pfp_shower_dedx', 'slc_pfp_track_score'
                  ]]
    temp["pfp_idx"] = temp.groupby(['run', 'subrun', 'event', 'slc_id']).transform("cumcount").add(1) - 1
    temp = temp.set_index(['run', 'subrun', 'event', 'slc_id']).reset_index()
    
    temp1 = temp.groupby(['run', 'subrun', 'event', 'slc_id'])['slc_pfp_shower_energy'].apply(list).reset_index()
    temp2 = temp.groupby(['run', 'subrun', 'event', 'slc_id'])['pfp_idx'].apply(list).reset_index()
    temp3 = temp.groupby(['run', 'subrun', 'event', 'slc_id'])['slc_pfp_shower_dir_x'].apply(list).reset_index()
    temp4 = temp.groupby(['run', 'subrun', 'event', 'slc_id'])['slc_pfp_shower_dir_y'].apply(list).reset_index()
    temp5 = temp.groupby(['run', 'subrun', 'event', 'slc_id'])['slc_pfp_shower_dir_z'].apply(list).reset_index()
    temp6 = temp.groupby(['run', 'subrun', 'event', 'slc_id'])['slc_pfp_shower_theta'].apply(list).reset_index()
    temp7 = temp.groupby(['run', 'subrun', 'event', 'slc_id'])['slc_pfp_shower_phi'].apply(list).reset_index()
    temp8 = temp.groupby(['run', 'subrun', 'event', 'slc_id'])['slc_pfp_shower_conv_gap'].apply(list).reset_index()
    temp9 = temp.groupby(['run', 'subrun', 'event', 'slc_id'])['slc_pfp_shower_dedx'].apply(list).reset_index()
    temp10 = temp.groupby(['run', 'subrun', 'event', 'slc_id'])['slc_pfp_track_score'].apply(list).reset_index()
    
    concat = pd.concat([temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10], axis = 1)
    
    concat = concat.loc[:, ~concat.columns.duplicated()]
    
    if len(concat) == 0:
        return df, df
    
    concat['n_pfp'] = concat.apply(lambda row: len(row['slc_pfp_shower_energy']), axis=1)
    
    #1 shower selection
    shw1 = concat[concat['n_pfp'] == 1]
    shw1 = shw1.drop(columns = ['slc_pfp_shower_energy', 'pfp_idx'
                                , 'slc_pfp_shower_dir_x', 'slc_pfp_shower_dir_y', 'slc_pfp_shower_dir_z'
                                ,'slc_pfp_shower_theta' ,'slc_pfp_shower_phi', 'slc_pfp_shower_conv_gap'
                                , 'slc_pfp_shower_dedx', 'slc_pfp_track_score'
                               ])
    shw1 = shw1.merge(temp, how='inner', on=['run','subrun','event','slc_id'])
    
    #2+ shower selection
    shw2more = concat[concat['n_pfp'] >= 2]
    temp = temp.drop(columns = ['slc_pfp_shower_energy', 'pfp_idx'
                                , 'slc_pfp_shower_dir_x', 'slc_pfp_shower_dir_y', 'slc_pfp_shower_dir_z'
                                ,'slc_pfp_shower_theta' ,'slc_pfp_shower_phi', 'slc_pfp_shower_conv_gap'
                                , 'slc_pfp_shower_dedx', 'slc_pfp_track_score'
                               ])
    temp = temp[~temp.duplicated(['run','subrun','event','slc_id'])]
    shw2more = shw2more.merge(temp, how='inner', on=['run','subrun','event','slc_id'])
    
    return shw1, shw2more
#------------------------------------------------------------------------------------------------------------------#   
def merge_df_prim(df, dfprim):
    dfprim = dfprim[['run','subrun','event','slc_id']]
    df = df.merge(dfprim, how='inner', on=['run','subrun','event','slc_id'])
    
    return df

#------------------------------------------------------------------------------------------------------------------#      
def ThetaFracFunc(theta, frac, m, c):
    return theta < ( frac + 1 ) * m + c

#------------------------------------------------------------------------------------------------------------------#   
def cutThetaFrac(df, m, c):
    df['thetaFrac'] = df.apply(lambda row: ThetaFracFunc(row['slc_pfp_shower_theta'], row['slc_opt0_frac'], m, c), axis = 1)
    
    df = df[df['thetaFrac'] == True]
    
    df = df.drop(['thetaFrac'], axis = 1)
    return df

#------------------------------------------------------------------------------------------------------------------#     
def vary_ThetaFrac(df_hnl, df_nu, df_cos, true_counts, start_counts):
    
    step = 1
    mStep = np.arange(5, 20 + step, step)
    step = 2
    cStep = np.arange(-20, 20 + step, step)
    
    p_arr, eff_arr, peff_arr = [], [], []
    m_arr, c_arr = [], []
    
    for m in mStep:
        for c in cStep:
            df_hnl_cut = cutThetaFrac(df_hnl, m, c) 
            df_nu_cut = cutThetaFrac(df_nu, m, c) 
            df_cos_cut = cutThetaFrac(df_cos, m, c) 
        
            purity, eff = calc_purity_eff(df_hnl_cut, df_nu_cut, df_cos_cut, true_counts, start_counts)
        
            p_arr.append(purity)
            eff_arr.append(eff)
            peff_arr.append(purity*eff)
            
            m_arr.append(m)
            c_arr.append(c)
    
    peffMax = max(peff_arr)
    bestIndexPE = peff_arr.index(peffMax)
    bestMPE = m_arr[bestIndexPE]
    bestCPE = c_arr[bestIndexPE]
    
    print("-------------------------------------")
    print("Best M PE = {0:3g}".format(bestMPE))
    print("Best C PE = {0:3g}".format(bestCPE))
    print("Purity = {0:3g}".format(p_arr[bestIndexPE]))
    print("Eff = {0:3g}".format(eff_arr[bestIndexPE]))

    pMax = max(p_arr)
    bestIndexP = p_arr.index(pMax)
    bestMP = m_arr[bestIndexP]
    bestCP = c_arr[bestIndexP]
    
    print("-------------------------------------")
    print("Best M PE = {0:3g}".format(bestMP))
    print("Best C PE = {0:3g}".format(bestCP))
    print("Purity = {0:3g}".format(p_arr[bestIndexP]))
    print("Eff = {0:3g}".format(eff_arr[bestIndexP]))
