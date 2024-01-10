import sys                                        
import os                                         
                                                  
os.nice(20)                                       
                                                  
import warnings                                   
warnings.filterwarnings("ignore")                 
                                                  
import argparse                                   
import uproot                                     
import h5py                                       
import math

# Local helper script                             
hnlDIR = os.environ['_']                          
sys.path.append(hnlDIR + '/pyscript/')            
                                                
from Plotting import *                            
from Dictionary import *
from HelperFunctions import *
from CutFunctions import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main(args):                                            

    #read in hdf5 to dataframe
    hnl_path = "./hdf5_files/hnl_test_precut_5k.h5"
    nu_path = "./hdf5_files/nu_test_precut_50k.h5"
    cosmics_path = "./hdf5_files/cosmics_test_precut_5k.h5"

    dfslc_hnl, dfsubrun_hnl, dfmevprtl_hnl, dfmct_hnl = hdf5_to_dataframe(hnl_path) 
    dfslc_nu, dfsubrun_nu, _, dfmct_nu = hdf5_to_dataframe(nu_path) 
    dfslc_cosmics, dfsubrun_cosmics, _, dfmct_cosmics = hdf5_to_dataframe(cosmics_path) 
    #----------------------------------------------------------------#

    #do scaling
    scale_pot_hnl, hnl_spill = calc_scaling_pot(dfsubrun_hnl, dfslc_hnl)
    scale_pot_nu, nu_spill = calc_scaling_pot(dfsubrun_nu, dfslc_nu)
    scale_pot_cosmics = calc_scaling_spill(dfsubrun_cosmics, dfslc_cosmics, hnl_spill, nu_spill)
    #----------------------------------------------------------------#
    
    #temp fix:: eventtype 0 in nu sample is unknown
    dfslc_nu['slc_true_event_type'][dfslc_nu['slc_true_event_type'] == 0] = -1
    
    #temp fix:: eventtype cosmics sample is -1 unknown --> change it to 9 
    dfslc_cosmics['slc_true_event_type'][dfslc_cosmics['slc_true_event_type'] == -1] = 9

    #dfslc_hnl['slc_true_event_type'][dfslc_hnl['slc_comp'] < 0.5] = -99
    #----------------------------------------------------------------#
    
    #convert opt0 unit from us to ns
    #dfslc_hnl["slc_opt0_time_corrected_Z_pandora"] = dfslc_hnl["slc_opt0_time_corrected_Z_pandora"] * 1000
    #dfslc_nu["slc_opt0_time_corrected_Z_pandora"] = dfslc_nu["slc_opt0_time_corrected_Z_pandora"] * 1000
    #dfslc_cosmics["slc_opt0_time_corrected_Z_pandora"] = dfslc_cosmics["slc_opt0_time_corrected_Z_pandora"] * 1000

    #----------------------------------------------------------------#
    #for efficiency calculation
    #true_counts, true_nonfv_counts = get_true_signal_in_all_spills(dfmct_hnl , scale_pot_hnl)
    #istart_counts, start_nonfv_counts = get_reco_signal_in_all_spills(dfslc_hnl, scale_pot_hnl)
   
    true_counts = 1019.6841605507992
    true_nonfv_counts = 434.88713106460574
    total_true_counts = 1454.571291615405
    start_counts = 957.0485396261358
    start_nonfv_counts = 362.1585664838355
    total_start_counts = 1319.2071061099714

    #true reco
    print("true signals = " + str(true_counts))
    print("true nonfv signals = " + str(true_nonfv_counts))
    print("total true signals = " +str(true_counts +true_nonfv_counts))
    
    #start reco signals
    print("start signals = " + str(start_counts))
    print("start nonfv signals = " + str(start_nonfv_counts))
    print("total start signals = " +str(start_counts + start_nonfv_counts))
    
    ##NO CUT---------------------------------------------------------------------
    #tag = '_nocut'

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'n_slc', 
    #                    tag,
    #                    xmin = 0, xmax = 50, xnbin = 50,
    #                    xtitle = 'nSlice'
    #                    )
    #
    ##PRE-SELECTION---------------------------------------------------------------------
    #tag = '_presel'

    #dfslc_hnl = cutPreSelection(dfslc_hnl)
    #dfslc_nu = cutPreSelection(dfslc_nu)
    #dfslc_cosmics = cutPreSelection(dfslc_cosmics)
    #
    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'n_slc', 
    #                    tag,
    #                    xmin = 0, xmax = 50, xnbin = 50,
    #                    xtitle = 'nSlice'
    #                    )
    #
    ##-------------------------------#
    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_crumbs_score', 
    #                    tag,
    #                    xmin = -1.5, xmax = 0.5, xnbin = 40,
    #                    xtitle = 'CRUMBS Score'
    #                    )

    #COSMICS-REJECTION---------------------------------------------------------------------
    #tag = '_cosrej'

    #dfslc_hnl = cutCosmics(dfslc_hnl)
    #dfslc_nu = cutCosmics(dfslc_nu)
    #dfslc_cosmics = cutCosmics(dfslc_cosmics)

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_crumbs_score', 
    #                    tag,
    #                    xmin = -1.5, xmax = 0.5, xnbin = 40,
    #                    xtitle = 'CRUMBS Score'
    #                    )
    #-------------------------------#
    #  
    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_track_end_x',  
    #                    tag,
    #                    xmin = -200, xmax = 200, xnbin = 100,
    #                    xtitle = 'PFParticle Track End X [cm]',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_track_end_y',  
    #                    tag,
    #                    xmin = -200, xmax = 200, xnbin = 100,
    #                    xtitle = 'PFParticle Track End Y [cm]',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_track_end_z',  
    #                    tag,
    #                    xmin = 0, xmax = 500, xnbin = 100,
    #                    xtitle = 'PFParticle Track End Z [cm]',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )
    ##CONTAINMENT---------------------------------------------------------------------
    tag = '_cont'

    dfslc_hnl = cutContainment(dfslc_hnl)
    dfslc_nu = cutContainment(dfslc_nu)
    dfslc_cosmics = cutContainment(dfslc_cosmics)

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_track_end_x',  
    #                    tag,
    #                    xmin = -200, xmax = 200, xnbin = 100,
    #                    xtitle = 'PFParticle Track End X [cm]',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_track_end_y',  
    #                    tag,
    #                    xmin = -200, xmax = 200, xnbin = 100,
    #                    xtitle = 'PFParticle Track End Y [cm]',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_track_end_z',  
    #                    tag,
    #                    xmin = 0, xmax = 500, xnbin = 100,
    #                    xtitle = 'PFParticle Track End Z [cm]',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )
    #
    #-------------------------------#
    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_opt0_score', 
    #                    tag,
    #                    xmin = 0, xmax = 40000, xnbin = 160,
    #                    xtitle = 'Opt0 Score'
    #                    )

    #OPT0 Score CUT---------------------------------------------------------------------
    tag = '_opt0score'

    dfslc_hnl = cutOpt0Score(dfslc_hnl)
    dfslc_nu = cutOpt0Score(dfslc_nu)
    dfslc_cosmics = cutOpt0Score(dfslc_cosmics)
   
    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_opt0_score', 
    #                    tag,
    #                    xmin = 0, xmax = 40000, xnbin = 160,
    #                    xtitle = 'Opt0 Score'
    #                    )

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts,
    #                    "slc_opt0_time_corrected_Z_pandora",
    #                    tag,
    #                    xmin = 1768, xmax = 1905, xnbin = 137,
    #                    xtitle = 'Opt0 Time Corrected Z Using Pandora Vertex [ns]',
    #                    ifPlotTime = True
    #                    )
    
    #TIME CUT---------------------------------------------------------------------
    tag = '_betweenBucket'

    dfslc_hnl = cutBetweenBucket(dfslc_hnl)
    dfslc_nu = cutBetweenBucket(dfslc_nu)
    dfslc_cosmics = cutBetweenBucket(dfslc_cosmics)
   
    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts,
    #                    "slc_opt0_time_corrected_Z_pandora",
    #                    tag,
    #                    xmin = 1768, xmax = 1905, xnbin = 137,
    #                    xtitle = 'Opt0 Time Corrected Z Using Pandora Vertex [ns]',
    #                    ifPlotTime = True
    #                    )

    #-------------------------------#
    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_n_razzled_muons', 
    #                    tag,
    #                    xmin = 0, xmax = 6, xnbin = 6,
    #                    xtitle = 'nRazzled Muons'
    #                    )

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_razzled_muon_score',  
    #                    tag,
    #                    xmin = 0, xmax = 1, xnbin = 50,
    #                    xtitle = 'PFParticle Razzled Muon Score',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )
    #MUON-REJECTION---------------------------------------------------------------------
    tag = '_muonrej'

    dfslc_hnl = cutMuon(dfslc_hnl)
    dfslc_nu = cutMuon(dfslc_nu)
    dfslc_cosmics = cutMuon(dfslc_cosmics)

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_n_razzled_muons', 
    #                    tag,
    #                    xmin = 0, xmax = 6, xnbin = 6,
    #                    xtitle = 'nRazzled Muons'
    #                    )

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_razzled_muon_score',  
    #                    tag,
    #                    xmin = 0, xmax = 1, xnbin = 50,
    #                    xtitle = 'PFParticle Razzled Muon Score',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )
    #-------------------------------#
    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_n_razzled_protons_thresh', 
    #                    tag,
    #                    xmin = 0, xmax = 6, xnbin = 6,
    #                    xtitle = 'nRazzled Protons'
    #                    )

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_razzled_proton_score',  
    #                    tag,
    #                    xmin = 0, xmax = 1, xnbin = 50,
    #                    xtitle = 'PFParticle Razzled Proton Score',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )
    #PROTON-REJECTION---------------------------------------------------------------------
    tag = '_protonrej'

    dfslc_hnl = cutProton(dfslc_hnl, ifScore =True)
    dfslc_nu = cutProton(dfslc_nu, ifScore =True)
    dfslc_cosmics = cutProton(dfslc_cosmics, ifScore =True)

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_n_razzled_protons_thresh', 
    #                    tag,
    #                    xmin = 0, xmax = 6, xnbin = 6,
    #                    xtitle = 'nRazzled Protons'
    #                    )

    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_razzled_proton_score',  
    #                    tag,
    #                    xmin = 0, xmax = 1, xnbin = 50,
    #                    xtitle = 'PFParticle Razzled Proton Score',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )
    ##-------------------------------#
    ##plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    ##                    true_counts, start_counts, 
    ##                    'slc_n_razzled_pions_thresh', 
    ##                    tag,
    ##                    xmin = 0, xmax = 6, xnbin = 6,
    ##                    xtitle = 'nRazzled Pions'
    ##                    )

    ##plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    ##                    true_counts, start_counts, 
    ##                    'slc_pfp_razzled_pion_score',  
    ##                    tag,
    ##                    xmin = 0, xmax = 1, xnbin = 50,
    ##                    xtitle = 'PFParticle Razzled Pion Score',
    ##                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    ##                    )

    ###PION-REJECTION---------------------------------------------------------------------
    #tag = '_pionrej'

    #dfslc_hnl = cutPion(dfslc_hnl,ifnPion = True, ifScore = True)
    #dfslc_nu = cutPion(dfslc_nu, ifnPion = True, ifScore = True)
    #dfslc_cosmics = cutPion(dfslc_cosmics, ifnPion=True, ifScore = True)
    #
    ##plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    ##                    true_counts, start_counts, 
    ##                    'slc_n_razzled_pions_thresh', 
    ##                    tag,
    ##                    xmin = 0, xmax = 6, xnbin = 6,
    ##                    xtitle = 'nRazzled Pions'
    ##                    )

    ##plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    ##                    true_counts, start_counts, 
    ##                    'slc_pfp_razzled_pion_score',  
    ##                    tag,
    ##                    xmin = 0, xmax = 1, xnbin = 50,
    ##                    xtitle = 'PFParticle Razzled Pion Score',
    ##                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    ##                    )
    ##-------------------------------#

    ##calculate angle w.r.t. beam
    #dfslc_hnl['slc_pfp_angle2Beam'] = dfslc_hnl.apply(lambda x: calc_angle2Beam(x['slc_pfp_shower_dir_x'], x['slc_pfp_shower_dir_y'], x['slc_pfp_shower_dir_z']), axis=1)
    #dfslc_nu['slc_pfp_angle2Beam'] = dfslc_nu.apply(lambda x: calc_angle2Beam(x['slc_pfp_shower_dir_x'], x['slc_pfp_shower_dir_y'], x['slc_pfp_shower_dir_z']), axis=1)
    #dfslc_cosmics['slc_pfp_angle2Beam'] = dfslc_cosmics.apply(lambda x: calc_angle2Beam(x['slc_pfp_shower_dir_x'], x['slc_pfp_shower_dir_y'], x['slc_pfp_shower_dir_z']), axis=1)

    plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
                        true_counts, start_counts, 
                        'slc_n_shws', 
                        tag,
                        xmin = 0, xmax = 6, xnbin = 6,
                        xtitle = 'nShowers'
                        )
    
    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_angle2Beam',  
    #                    tag,
    #                    xmin = 0, xmax = 180, xnbin = 180,
    #                    xtitle = 'PFParticle Shower Angle To Beam [degrees]',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )

    plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
                        true_counts, start_counts, 
                        'slc_pfp_razzled_electron_score',  
                        tag,
                        xmin = 0, xmax = 1, xnbin = 50,
                        xtitle = 'PFParticle Razzled Electron Score',
                        ytitle = 'PFParticles (1x10$^{21}$ POT)'
                        )

    plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
                        true_counts, start_counts, 
                        'slc_n_razzled_electrons', 
                        tag,
                        xmin = 0, xmax = 6, xnbin = 6,
                        xtitle = 'nRazzled Electrons'
                        )

    plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
                        true_counts, start_counts, 
                        'slc_pfp_razzled_photon_score',  
                        tag,
                        xmin = 0, xmax = 1, xnbin = 50,
                        xtitle = 'PFParticle Razzled Photon Score',
                        ytitle = 'PFParticles (1x10$^{21}$ POT)'
                        )

    plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
                        true_counts, start_counts, 
                        'slc_n_razzled_photons', 
                        tag,
                        xmin = 0, xmax = 6, xnbin = 6,
                        xtitle = 'nRazzled Photons'
                        )

    ##SHOWER-REQUIREMENT---------------------------------------------------------------------
    tag = '_shwreq'

    dfslc_hnl = cutShower(dfslc_hnl, ifDir = False)
    dfslc_nu = cutShower(dfslc_nu, ifDir = False)
    dfslc_cosmics = cutShower(dfslc_cosmics, ifDir = False)

    plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
                        true_counts, start_counts, 
                        'slc_n_shws', 
                        tag,
                        xmin = 0, xmax = 6, xnbin = 6,
                        xtitle = 'nShowers'
                        )
    
    #plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
    #                    true_counts, start_counts, 
    #                    'slc_pfp_angle2Beam',  
    #                    tag,
    #                    xmin = 0, xmax = 180, xnbin = 180,
    #                    xtitle = 'PFParticle Shower Angle To Beam [degrees]',
    #                    ytitle = 'PFParticles (1x10$^{21}$ POT)'
    #                    )

    plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
                        true_counts, start_counts, 
                        'slc_pfp_razzled_electron_score',  
                        tag,
                        xmin = 0, xmax = 1, xnbin = 50,
                        xtitle = 'PFParticle Razzled Electron Score',
                        ytitle = 'PFParticles (1x10$^{21}$ POT)'
                        )

    plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
                        true_counts, start_counts, 
                        'slc_n_razzled_electrons', 
                        tag,
                        xmin = 0, xmax = 6, xnbin = 6,
                        xtitle = 'nRazzled Electrons'
                        )

    plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
                        true_counts, start_counts, 
                        'slc_pfp_razzled_photon_score',  
                        tag,
                        xmin = 0, xmax = 1, xnbin = 50,
                        xtitle = 'PFParticle Razzled Photon Score',
                        ytitle = 'PFParticles (1x10$^{21}$ POT)'
                        )

    plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
                        true_counts, start_counts, 
                        'slc_n_razzled_photons', 
                        tag,
                        xmin = 0, xmax = 6, xnbin = 6,
                        xtitle = 'nRazzled Photons'
                        )


if __name__ == '__main__':                                 
                                                           
    parser = argparse.ArgumentParser()                     
    args = parser.parse_args()                             
    main(args)                                             
