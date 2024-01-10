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
    hnl_path = "./hdf5_files/hnl_1k.h5"
    
    dfslc_hnl, _, _, _ = hdf5_to_dataframe(hnl_path) 

    whenHNL = dfslc_hnl['slc_true_event_type'] == 0
    hnl_true_T = dfslc_hnl['slc_true_vtx_t'][whenHNL] 
   
    c_cm_per_ns = 29.9792458 
    hnl_true_z = dfslc_hnl['slc_true_vtx_z'][whenHNL] 
    hnl_true_Tcorr = hnl_true_T - hnl_true_z/c_cm_per_ns

    whenCOS = dfslc_hnl['slc_true_event_type'] == 9
    cos_true_T = dfslc_hnl['slc_true_vtx_t'][whenCOS]

    #========================
    xmin = 360
    xmax = 660
    xnbin = 300
    bins = np.arange(xmin, xmax + 2*(xmax-xmin)/xnbin, (xmax-xmin)/xnbin),

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize = (16, 12))
    
    plot_1dhist(
                hnl_true_T, 
                ax1,
                xmin, xmax, xnbin,
                xlimmin = xmin, xlimmax = xmax,
		ifnorm = True,
		label = 'hnl', color = 'orange', iflabelbox = True
                )
    
    plot_1dhist(
                hnl_true_Tcorr, 
                ax2,
                xmin, xmax, xnbin,
                xlimmin = xmin, xlimmax = xmax,
		ifnorm = True,
		label = 'hnl corr', color = 'green', iflabelbox = True
                )
    
    plot_1dhist(
                cos_true_T, 
                ax3,
                xmin, xmax, xnbin,
                xlimmin = xmin, xlimmax = xmax,
		ifnorm = True,
		label = 'cosmics', color = 'blue', iflabelbox = True
                )

    plt.savefig("./plot_files/test_true_time.png", dpi = 200)


if __name__ == '__main__':                                 
                                                           
    parser = argparse.ArgumentParser()                     
    args = parser.parse_args()                             
    main(args)                                             
