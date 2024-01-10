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
    hnl_path = "./hdf5_files/hnl_5k.h5"
    nu_path = "./hdf5_files/nu_50k.h5"
    cosmics_path = "./hdf5_files/cosmics_5k.h5"

    dfslc_hnl, dfsubrun_hnl, dfmevprtl_hnl, dfmct_hnl = hdf5_to_dataframe(hnl_path) 
    dfslc_nu, dfsubrun_nu, _, dfmct_nu = hdf5_to_dataframe(nu_path) 
    dfslc_cosmics, dfsubrun_cosmics, _, dfmct_cosmics = hdf5_to_dataframe(cosmics_path) 
    #----------------------------------------------------------------#

    #do scaling
    scale_pot_hnl, hnl_spill = calc_scaling_pot(dfsubrun_hnl, dfslc_hnl)
    scale_pot_nu, nu_spill = calc_scaling_pot(dfsubrun_nu, dfslc_nu)
    print(hnl_spill + nu_spill)
    scale_pot_cosmics = calc_scaling_spill(dfsubrun_cosmics, dfslc_cosmics, hnl_spill, nu_spill)
    #----------------------------------------------------------------#
    
    #temp fix:: eventtype 0 in nu sample is unknown
    dfslc_nu['slc_true_event_type'][dfslc_nu['slc_true_event_type'] == 0] = -1
    
    #temp fix:: eventtype cosmics sample is -1 unknown --> change it to 9 
    dfslc_cosmics['slc_true_event_type'][dfslc_cosmics['slc_true_event_type'] == -1] = 9

#    #dfslc_hnl['slc_true_event_type'][dfslc_hnl['slc_comp'] < 0.5] = -99
#    #----------------------------------------------------------------#
#    #for efficiency calculation
#    true_counts, true_nonfv_counts = get_true_signal_in_all_spills(dfmct_hnl , scale_pot_hnl)
#    start_counts, start_nonfv_counts = get_reco_signal_in_all_spills(dfslc_hnl, scale_pot_hnl)
#    
#    #true reco
#    print("true signals = " + str(true_counts))
#    print("true nonfv signals = " + str(true_nonfv_counts))
#    print("total true signals = " +str(true_counts +true_nonfv_counts))
#    
#    #start reco signals
#    print("start signals = " + str(start_counts))
#    print("start nonfv signals = " + str(start_nonfv_counts))
#    print("total start signals = " +str(start_counts + start_nonfv_counts))
#    
#    #NO CUT---------------------------------------------------------------------
#    tag = '_nocut'
#
#    plot_slc_var(dfslc_hnl, dfslc_nu, dfslc_cosmics,
#                        true_counts, start_counts, 
#                        'n_slc', 
#                        tag,
#                        xmin = 0, xmax = 50, xnbin = 50,
#                        xtitle = 'nSlice'
#                        )
    

if __name__ == '__main__':                                 
                                                           
    parser = argparse.ArgumentParser()                     
    args = parser.parse_args()                             
    main(args)                                             
