import os
import sys
os.nice(20)

import math
from math import floor, log10
import pandas as pd
from matplotlib.patches import Rectangle

# Local helper script                             
hnlDIR = os.environ['_']                          
sys.path.append(hnlDIR + '/pyscript/')            
                                                
from Plotting import *                            
from Dictionary import *

fontsize = 12

#------------------------------------------------------------------------------------------------------------------#
# Define function for string formatting of scientific notation
def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if exponent is None:
        exponent = int(floor(log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

    return r"${0:.{2}f}\times10^{{{1:d}}}$".format(coeff, exponent, precision)

#------------------------------------------------------------------------------------------------------------------#
def plot_purity_eff(p_arr, eff_arr, peff_arr, cutStep,  purity_start, eff_start, loc):
    
    peffMax = max(peff_arr)
    bestIndexPE = peff_arr.index(peffMax)
    bestScorePE = cutStep[bestIndexPE]
    
    print("-------------------------------------")
    print("Best Cut Score PE = {0:3g}".format(bestScorePE))
    print("Purity = {0:3g}".format(p_arr[bestIndexPE]))
    print("Eff = {0:3g}".format(eff_arr[bestIndexPE]))
    print("Purity - Start Purity = {0:3g}".format(p_arr[bestIndexPE] - purity_start))
    print("Eff - Start Eff = {0:3g}".format(eff_arr[bestIndexPE] - eff_start))
    
    pMax = max(p_arr)
    bestIndexP = p_arr.index(pMax)
    bestScoreP = cutStep[bestIndexP]
    
    print("-------------------------------------")
    print("Best Cut Score P = {0:3g}".format(bestScoreP))
    print("Purity = {0:3g}".format(p_arr[bestIndexP]))
    print("Eff = {0:3g}".format(eff_arr[bestIndexP]))
    print("Purity - Start Purity = {0:3g}".format(p_arr[bestIndexP] - purity_start))
    print("Eff - Start Eff = {0:3g}".format(eff_arr[bestIndexP] - eff_start))
    
    fig, ax1 = plt.subplots(1, 1, figsize=(6, 4), sharex = True)
    ax2 = ax1.twinx()

    lns1 = ax1.plot(cutStep, eff_arr, c='g')
    lns2 = ax2.plot(cutStep, p_arr, c='orange')
    lns3 = ax1.axvline(x = bestScorePE, color = 'r', label = "Best Purity*Select Eff = {0:4g}".format(bestScorePE))
    lns4 = ax1.axvline(x = bestScoreP, color = 'purple', label = "Best Purity Cut = {0:4g}".format(bestScoreP))

    #ax1.set_xlabel("Cut Score")
    ax1.set_xlabel("Cut Score")
    ax1.set_ylabel("Selection Efficiency [%]", c='g')
    ax1.set_ylim(-5, 115)
    ax2.set_ylabel("Purity [%]", c = 'orange')
    ax2.set_ylim(-5, 115)
    
    handles, labels = ax1.get_legend_handles_labels()

    labels.insert(1,"Purity = {0:.3g}%".format(p_arr[bestIndexPE]))
    labels.insert(2,"Select Eff = {0:.3g}%".format(eff_arr[bestIndexPE]))
    labels.insert(3,"Eff*Purity = {0:.3g}%".format(p_arr[bestIndexPE]*eff_arr[bestIndexPE]))
                  
    labels.insert(5,"Purity = {0:.3g}%".format(p_arr[bestIndexP]))
    labels.insert(6,"Select Eff = {0:.3g}%".format(eff_arr[bestIndexP]))
    labels.insert(7,"Eff*Purity = {0:.3g}%".format(p_arr[bestIndexP]*eff_arr[bestIndexP]))

    extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
    handles.insert(1, extra)
    handles.insert(2, extra)
    handles.insert(3, extra)
    handles.insert(5, extra)
    handles.insert(6, extra)
    handles.insert(7, extra)

    ax1.legend(handles, labels, bbox_to_anchor=(0.65, 0.23), fontsize=fontsize - 4, fancybox=False, ncol = 1)
    ax1.legend(handles, labels, loc=loc, fontsize=fontsize - 4, fancybox=False, ncol = 1)

#------------------------------------------------------------------------------------------------------------------#
def getUfromScaleFactor(inputU, scaleFactor):
    return np.sqrt(scaleFactor)*inputU
#------------------------------------------------------------------------------------------------------------------#
def make_interval(lb_val, ub_val):

    lb_arr = []
    ub_arr = []   
    
    #construct the lower/upper bound array
    for idx in range(0,81):
        lb_arr.append(lb_val+19*idx)
        ub_arr.append(ub_val+19*idx)
    return lb_arr, ub_arr

def checkInterval(x, lb_arr, ub_arr):

    isIn = False

    for lb, ub in zip(lb_arr, ub_arr):
        if lb <= x <= ub:
            isIn = True
            break

    return isIn
#------------------------------------------------------------------------------------------------------------------#
def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    angle = math.degrees(angle)

    return angle

def calc_angle2Beam(x, y, z):

    v1 = [0, 0, 1]
    v2 = [x, y, z]
    angle = angle_between(v1, v2)
    
    return angle

#------------------------------------------------------------------------------------------------------------------#
def check_dataframe(df):
    print(df.shape)
    print(df.columns)
    print(df.head(20))
    print("--------------------------------")

#------------------------------------------------------------------------------------------------------------------#
def calc_scaling_pot(df, dfslc, ifScale = 1):
    target_pot = 1*10**21
    pot_per_spill = 5e12

    sample_pot = df['pot'].sum()
    sample_spill = df['ngenevts'].sum()
    
    if ifScale != 1:
        sample_pot = sample_pot / ifScale
    
    scale_pot = target_pot / sample_pot
    
    target_spill = target_pot * sample_spill / sample_pot

    dfslc['scale_pot'] = scale_pot
    
    print('-----------------------------------------------')
    print('sample pot = {:2f}'.format(sample_pot))
    print('sample spill = {:2f}'.format(sample_spill))
    print('target spill = {:2f}'.format(target_spill))
    print('scale pot factor = {:2f}'.format(scale_pot))
    print('-----------------------------------------------')

    return scale_pot, target_spill

#------------------------------------------------------------------------------------------------------------------#
def calc_scaling_spill(df, dfslc, hnl_spill, nu_spill):

    target_pot = 1*10**21
    pot_per_spill = 5e12

    target_spill = target_pot / pot_per_spill 
    target_intime_spill = target_spill - hnl_spill - nu_spill

    sample_spill = df['ngenevts'].sum()

    scale_pot = target_intime_spill / sample_spill
    
    dfslc['scale_pot'] = scale_pot

    print('-----------------------------------------------')
    print('target total spill = ' + str(target_spill))
    print('hnl + nu spill = ' + str(hnl_spill + nu_spill))
    print('target intime spill = ' + str(target_intime_spill))
    print('scale pot factor = ' + str(scale_pot))
    print('-----------------------------------------------')

    return scale_pot 

#------------------------------------------------------------------------------------------------------------------#
def get_true_signal_in_all_spills(df, scale_pot):
    
    nSig = len(df['nu_event_type'][df['nu_event_type']==0]) * scale_pot
    nNonFV = len(df['nu_event_type'][df['nu_event_type']==1]) * scale_pot
    
    return nSig, nNonFV

#------------------------------------------------------------------------------------------------------------------#
def get_reco_signal_in_all_spills(df, scale_pot):
    
    #dataframe are at pfp level --> need to count slice
    df = df.drop_duplicates(subset=["run","subrun","event", "slc_id"])
    
    #keep only relevant columns
    df = df[['run', 'subrun', 'event', 'slc_id', 'slc_comp', 'slc_true_event_type']]
    
    #count
    nSig = len(df['slc_true_event_type'][ (df['slc_true_event_type']==0) & (df['slc_comp'] > 0.5)]) * scale_pot
    nNonFV = len(df['slc_true_event_type'][ (df['slc_true_event_type']==1) & (df['slc_comp'] > 0.5)]) * scale_pot

    return nSig, nNonFV

#------------------------------------------------------------------------------------------------------------------#
def get_slice(df, scale_pot, ifSignal):
   
    #dataframe are at pfp level --> need to count slice
    df = df.drop_duplicates(subset=["run","subrun","event", "slc_id"])
    
    #keep only relevant columns
    df = df[['run', 'subrun', 'event', 'slc_id', 'slc_comp']]
   
    if ifSignal:
        return len(df[df['slc_comp'] > 0.5]) * scale_pot
    else:
        return len(df['slc_comp']) * scale_pot

#------------------------------------------------------------------------------------------------------------------#
def split_into_event_type(df, var_name):
    
    pltdf, potdf = [], []
    
    for t in event_type:

        #event type condition
        when = df['slc_true_event_type'] == t
        
        #plot variable - can be pfp, do not use for counting
        dt = df[var_name][when]
        pltdf.append(dt)
    
        #pot varible - make it the same shape as plotting variable
        dp = df['scale_pot'][when]
        potdf.append(dp)
        
    return pltdf, potdf

#------------------------------------------------------------------------------------------------------------------#
def count_slice(dfhnl, dfnu, dfcosmics, 
                true_counts, start_counts
               ):

    scale_pot_hnl = dfhnl['scale_pot'].unique()
    scale_pot_nu = dfnu['scale_pot'].unique()
    scale_pot_cosmics = dfcosmics['scale_pot'].unique()

    if len(scale_pot_hnl) == 0:
        scale_pot_hnl = 0
    else:
        scale_pot_hnl = scale_pot_hnl[0]  

    if len(scale_pot_cosmics) == 0:
        scale_pot_cosmics = 0
    else:
        scale_pot_cosmics = scale_pot_cosmics[0]

    slc_count_hnl, slc_count_nu, slc_count_cosmics = [], [], []
    n_hnl, n_nu, n_cosmics = 0, 0, 0

    for t in event_type:
        
        #event type condition
        whenhnl = dfhnl['slc_true_event_type'] == t
        whennu = dfnu['slc_true_event_type'] == t
        whencosmics = dfcosmics['slc_true_event_type'] == t
        
        if t == 0 or t == 1:
            n_hnl = get_slice(dfhnl[whenhnl], scale_pot_hnl, True)
        else:
            n_hnl = get_slice(dfhnl[whenhnl], scale_pot_hnl, False)
            n_cosmics = get_slice(dfcosmics[whencosmics], scale_pot_cosmics, False)
            
            #Neutrino sample has 2 POT so gotta be careful
            n_nu = 0
            if len(scale_pot_nu) == 0:
                n_nu += get_slice(dfnu[whennu], 0, False)
            else:
                for pot in scale_pot_nu:
                    whennuPOT = dfnu['scale_pot'] == pot  
                    n_nu += get_slice(dfnu[whennu & whennuPOT], pot, False)
                    #print("type {}, pot {}, non-scaled n_nu {}".format(event_label_dict[t], pot, get_slice(dfnu[whennu & whennuPOT], pot, False)/pot))
            
        slc_count_hnl.append(n_hnl)
        slc_count_nu.append(n_nu)
        slc_count_cosmics.append(n_cosmics)

    slc_count = [i + j + k for (i,j,k) in zip(slc_count_hnl, slc_count_nu, slc_count_cosmics)]

    total_hnl = slc_count_hnl[0] + slc_count_hnl[1]
    
    if sum(slc_count) > 0:
        purity = total_hnl / sum(slc_count) * 100
    else:
        purity = 0 
    
    eff = total_hnl / true_counts * 100
    
    select_eff = total_hnl / start_counts * 100

    total_bkg = sum(slc_count) - total_hnl
    
    update_label = [i + " (" + '{:,}'.format(round(j)) + ")" for (i,j) in zip(event_label, slc_count)]

    #print('nSig = {0:.6g}, nBkg = {1:.6g}, nSlc = {2:.6g}'.format(total_hnl, total_bkg, sum(slc_count)))
    #print('purity = {0:.3g}'.format(purity))
    #print('eff = {0:.3g}'.format(eff))
    #print('select eff = {0:.3g}'.format(select_eff))

    return update_label, purity, eff, select_eff

#------------------------------------------------------------------------------------------------------------------#
def calc_purity_eff(dfhnl, dfnu, dfcosmics,
                true_counts, start_counts 
                ):

    #keep only relevant columns
    dfhnl = dfhnl[['run', 'subrun', 'event', 'slc_id', 'slc_comp', 'slc_true_event_type', 'scale_pot']]
    dfnu = dfnu[['run', 'subrun', 'event', 'slc_id', 'slc_comp', 'slc_true_event_type', 'scale_pot']]
    dfcosmics = dfcosmics[['run', 'subrun', 'event', 'slc_id', 'slc_comp', 'slc_true_event_type', 'scale_pot']]

    #count slices and add to labels
    _, purity, _, select_eff = count_slice(dfhnl, dfnu, dfcosmics, true_counts, start_counts)
   
    purity = round(purity,1)
    select_eff = round(select_eff,1)
    return purity, select_eff

    
#------------------------------------------------------------------------------------------------------------------#
def plot_slc_var(dfhnl, dfnu, dfcosmics,
                true_counts, start_counts, 
                var_name, 
                xmin, xmax, xnbin,
                xtitle,
                ytitle =  "Slices (1x10$^{21}$ POT)",
                ifPlotTime = False,
                ifAddLegend = False,
                addLegend = "test", 
                LegendLoc = "best"
                ):

    #keep only relevant columns
    dfhnl = dfhnl[['run', 'subrun', 'event', 'slc_id', 'slc_comp', 'slc_true_event_type', 'scale_pot', var_name]]
    dfnu = dfnu[['run', 'subrun', 'event', 'slc_id', 'slc_comp', 'slc_true_event_type', 'scale_pot', var_name]]
    dfcosmics = dfcosmics[['run', 'subrun', 'event', 'slc_id', 'slc_comp', 'slc_true_event_type', 'scale_pot', var_name]]
    
    #count slices and add to labels
    update_label, purity, eff, select_eff = count_slice(dfhnl, dfnu, dfcosmics, true_counts, start_counts)
   
    # split the df into event type for plotting
    pltdf_hnl, potdf_hnl = split_into_event_type(dfhnl, var_name) 
    pltdf_nu, potdf_nu = split_into_event_type(dfnu, var_name)
    pltdf_cosmics, potdf_cosmics = split_into_event_type(dfcosmics, var_name)

    #bin hnl/nu separately for different pot scaling
    hist_hnl, bins_hnl, _ = plt.hist(pltdf_hnl, bins = xnbin, range =(xmin, xmax), weights = potdf_hnl)
    hist_nu, _, _ = plt.hist(pltdf_nu, bins = xnbin, range =(xmin, xmax), weights = potdf_nu)
    hist_cosmics, _, _ = plt.hist(pltdf_cosmics, bins = xnbin, range =(xmin, xmax), weights = potdf_cosmics)

    plt.clf()
    
    #add them back together again
    hist = np.add(hist_hnl, hist_nu)
    hist = np.add(hist, hist_cosmics)

    del hist_hnl, hist_nu, hist_cosmics
    
    #fake bottom for stacking histogram
    bottom = np.zeros(xnbin)

    fig, ax = plt.subplots(1, 1, figsize = (6, 4))

    #Plot Time Bucket----------------------------------------------------#
    if ifPlotTime:
        fig, ax = plt.subplots(1, 1, figsize = (12, 4))
        
        lb_arr, ub_arr = make_interval(375, 380)
        #lb_arr, ub_arr = make_interval(377, 380)
        
        for i in range(74, 81):
            ax.axvline(lb_arr[i], color = 'indigo', alpha = 0.5, lw = 1, linestyle = 'dotted')
            ax.axvline(ub_arr[i], color = 'indigo', alpha = 0.5, lw = 1, linestyle = 'dotted')
            ax.axvspan(lb_arr[i], ub_arr[i], alpha=0.2, color='indigo')
    #----------------------------------------------------#

    for i in reversed(range(0, len(event_type))):
        plot_bar(
                 bins_hnl[:-1], hist[i],
                 ax,
                 width=np.diff(bins_hnl),
                 xlimmin = xmin, xlimmax = xmax,
                 bottom = bottom,
                 label = update_label[i],
                 color = event_col[i],
                 xtitle = xtitle,
                 ytitle = ytitle,
                 fontsize = fontsize,
                 edgecolor='white',
                 linewidth=0.2
                 )
        bottom += hist[i]

    ax.yaxis.set_major_locator(MaxNLocator(5, prune='lower'))
    
    handles, labels = plt.gca().get_legend_handles_labels()
    handles = handles[::-1]
    labels = labels[::-1]

    labels.append("Purity = {0:.3g}%".format(purity))
    labels.append("Eff = {0:.3g}%".format(eff))
    labels.append("Select Eff = {0:.3g}%".format(select_eff))

    extra = Rectangle((0, 0), 0.1, 0.1, fc="w", fill=False, edgecolor='none', linewidth=0)
    handles.append(extra)
    handles.append(extra)
    handles.append(extra)
    
    if ifAddLegend == True:
        #labels.insert(0, addLegend)
        #handles.insert(0, extra)
        num =  labels[0][labels[0].find("(")+1:labels[0].find(")")]
        num_str = num.split(",")[0] + num.split(",")[1]
        num_float = float(num_str)
        labels[0] = addLegend + '\n ({:,})'.format(round(num_float))

    ax.legend(handles, labels, loc=LegendLoc, fontsize=fontsize - 4, fancybox=False, ncol = 1)

    #ax.set_yscale('log')

    fig.tight_layout()
    
    #plt.savefig("./plot_files/" + var_name + plot_tag + ".png", dpi = 200)
    return hist, bins_hnl
