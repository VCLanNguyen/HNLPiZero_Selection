import os
import sys
os.nice(20)

#import ROOT
import array

# Local helper script                             
hnlDIR = os.environ['_']                          
sys.path.append(hnlDIR + '/pyscript/')            
                                                
from Plotting import *                            
from Dictionary import *

fontsize = 12
#------------------------------------------------------------------------------------------------------------------#
g4_list = ['slc_geant4_multisim_reinteractions']
g4_name = ['G4 Re-Interaction']

unisim_list = [ 'slc_xsec_unisim_DecayAngMEC', 
               'slc_xsec_unisim_ThetaDelta2NRad',
               'slc_xsec_unisim_Theta_Delta2Npi', 
               'slc_xsec_unisim_VecFFCCQEshape',]

multisigma_list = [
        'slc_xsec_multisigma_CoulombCCQE',
       'slc_xsec_multisigma_NonRESBGvbarnCC1pi',
       'slc_xsec_multisigma_NonRESBGvbarnCC2pi',
       'slc_xsec_multisigma_NonRESBGvbarnNC1pi',
       'slc_xsec_multisigma_NonRESBGvbarnNC2pi',
       'slc_xsec_multisigma_NonRESBGvbarpCC1pi',
       'slc_xsec_multisigma_NonRESBGvbarpCC2pi',
       'slc_xsec_multisigma_NonRESBGvbarpNC1pi',
       'slc_xsec_multisigma_NonRESBGvbarpNC2pi',
       'slc_xsec_multisigma_NonRESBGvnCC1pi',
       'slc_xsec_multisigma_NonRESBGvnCC2pi',
       'slc_xsec_multisigma_NonRESBGvnNC1pi',
       'slc_xsec_multisigma_NonRESBGvnNC2pi',
       'slc_xsec_multisigma_NonRESBGvpCC1pi',
       'slc_xsec_multisigma_NonRESBGvpCC2pi',
       'slc_xsec_multisigma_NonRESBGvpNC1pi',
       'slc_xsec_multisigma_NonRESBGvpNC2pi', 
       'slc_xsec_multisigma_NormCCMEC',
       'slc_xsec_multisigma_NormNCMEC', 
       'slc_xsec_multisigma_RDecBR1eta',
       'slc_xsec_multisigma_RDecBR1gamma', 
       'slc_xsec_multisigma_RPA_CCQE',
       'slc_xsec_multisigma_NormNCCOH', 
       'slc_xsec_multisigma_NormCCCOH'
]

multisim_list = [
        "slc_xsec_multisim_ZExpA",
        "slc_xsec_multisim_NCEL",
        "slc_xsec_multisim_CCRES",
        "slc_xsec_multisim_NCRES",
        "slc_xsec_multisim_DISBY",
        "slc_xsec_multisim_FSI_pi",
        "slc_xsec_multisim_FSI_N"
]

flux_list = [
        "slc_flux_weight_expskin",
        "slc_flux_weight_horncurrent",
        "slc_flux_weight_kminus",
        "slc_flux_weight_kplus",
        "slc_flux_weight_kzero",
        "slc_flux_weight_nucleoninexsec",
        "slc_flux_weight_nucleonqexsec",
        "slc_flux_weight_nucleontotxsec",
        "slc_flux_weight_piminus",
        "slc_flux_weight_pioninexsex",
        "slc_flux_weight_pionqexsec",
        "slc_flux_weight_piontotxsec",
        "slc_flux_weight_piplus"    
]

flux_name = ['Exposure Skin Weight'
                ,'Horn Current Weight'
                ,'Kaon Minus Weight'
                ,'Kaon Plus Weight'
                ,'Neutral Kaon Weight'
                ,'Nucleon Ineslastic Cross Section Weight'
                ,'Nucleon Quasi-Elastic Cross Section Weight'
                ,'Nucleon Total Cross Section Weight'
                ,'Pion Minus Weight'
                ,'Pion Inelastic Cross Section Weight'
                ,'Pion Quasi-Elastic Cross Section Weight'
                ,'Pion Total Cross Section Weight'
                ,'Pion Plus Weight'
                ]

#------------------------------------------------------------------------------------------------------------------#
nu_col = col_dict['Teal']
hnl_col = col_dict['Flamingo']
cv_col = col_dict['Coral']


total_col = col_dict['Purple']

stats_col = col_dict['Aqua']
flx_col = col_dict['MintGreen']

xsec_col = col_dict['Peach']
g4_col = col_dict['RosyBrown4']

#------------------------------------------------------------------------------------------------------------------#
xmin = 0
xmax = 19
xnbin = 19

hnl_ymin = 0
hnl_ymax = 2500
hnl_ymax2 = 0.2

rockbox_ymin = 0
rockbox_ymax = 50
rockbox_ymax2 = 2

ncpi0_ymin = 0
ncpi0_ymax = 400
ncpi0_ymax2 = 2

nu_ymin = 0
nu_ymax = 3000
nu_ymax2 = 2

cos_ymin = 0
cos_ymax = 8
cos_ymax2 = 2

bins = np.arange(xmin, xmax+(xmax-xmin)/xnbin, (xmax-xmin)/xnbin)
bins_mid = np.convolve(bins, [0.5, 0.5], "valid")

#------------------------------------------------------------------------------------------------------------------#
def fill_nan_plz(which_dict):
    for k,v in zip(which_dict.keys(), which_dict.values()):
        v = np.nan_to_num(v)
        which_dict[k] = v
#------------------------------------------------------------------------------------------------------------------#
def new_19_by_19_cov():
    w, h = 19, 19
    cov = [[0] * w for i in range(h)]
    cov = np.array(cov)
    
    return cov
#------------------------------------------------------------------------------------------------------------------#
def plot_hatchy_hatch(which_dict, label, which_type, error_type):
    
    ymin, ymax = 0,0
    #-----------------------------------------------------------------#
    if which_type == 'hnl':
        ymin = hnl_ymin
        ymax = hnl_ymax
    elif which_type == 'rockbox':
        ymin = rockbox_ymin
        ymax = rockbox_ymax
    elif which_type == 'ncpi0':
        ymin = ncpi0_ymin
        ymax = ncpi0_ymax
    elif which_type == 'cos':
        ymin = cos_ymin
        ymax = cos_ymax
        
    if error_type == 'stat_err':
        error_label = 'Statistics'
        col = stats_col
    elif error_type == 'flx_err':
        error_label = 'Flux'
        col = flx_col    
    elif error_type == 'xsec_err':
        error_label = 'Cross Section'
        col = xsec_col
    elif error_type == 'g4_err':
        error_label = 'Re-Interaction'
        col = g4_col
    #-----------------------------------------------------------------#
    fig, (ax1) = plt.subplots(1,1, figsize = (6,4))

    #nStat NoScale
    ax1.step(bins, which_dict['cv_plot'] 
         , color = cv_col
         , label = label
        )

    #-----------------------------------------------------------------#
    #stats

    bottom = which_dict['cv'] - which_dict[error_type]
    peak = which_dict['cv'] + which_dict[error_type]
    height = peak - bottom

    plt. rcParams["hatch.color"] = col
    ax1.bar(
        x = bins[:-1]
        , height= height
        , width = np.diff(bins)
        , bottom = bottom
        , align = 'edge'
        , hatch='///'
        , fc=col
        , alpha = 0.2 
        , label = error_label
            )
    #-----------------------------------------------------------------#

    ax1.legend(loc = 'upper left',fontsize = 14)
    plot_tick(ax1, 16)
    plot_title(ax1, "", 'Opt0 Time Corrected Z % 18.936 [ns]',  "Slices (No Scaling)", 16)

    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)

    #-----------------------------------------------------------------#
    fig.tight_layout()
    
#------------------------------------------------------------------------------------------------------------------#    
def get_cov_corr_matrix(cv, universe):

    print(cv.shape)
    print(universe.shape)
    
    #construct covariance matrix
    cov = np.cov(np.subtract(cv[np.newaxis,:], universe).transpose())

    #error = sqrt(diag)
    err = np.sqrt(np.diagonal(cov))

    #construct correlation matrix
    corr = cov/np.outer(err, err)  
    
    return cov, corr

#------------------------------------------------------------------------------------------------------------------#   
def make_df_multisim(df, name):
    
    #explode array into columns
    df_flxw = pd.DataFrame(df[name].tolist(),index=df.index).add_prefix(name+'_')

    #merge 2 arrays together
    df = pd.concat([df, df_flxw], axis = 1)
    
    del df_flxw
    return df

#------------------------------------------------------------------------------------------------------------------# 
def plot_and_save_universe(df, which_dict, name, good_name, len_univ, xmin, xmax, xnbin , which_type, savePath):
    universes = []

    ymin, ymax = 0,0
    #-----------------------------------------------------------------#
    if which_type == 'hnl':
        ymin = hnl_ymin
        ymax = hnl_ymax
    elif which_type == 'rockbox':
        ymin = rockbox_ymin
        ymax = rockbox_ymax
    elif which_type == 'ncpi0':
        ymin = ncpi0_ymin
        ymax = ncpi0_ymax
    #-----------------------------------------------------------------#
    fig, ax = plt.subplots(1,1, figsize = (6,4))

    #-----------------------------------------------------------------#
    #Universe
    pltdf = df['mod_t']
    
    for idx in range(0, len_univ):
        
        weights = df[name + '_{}'.format(idx)]
        label = ''
        if idx == 0:
            label = good_name + "\nSystematic Universes"
        
        univ, _, _ = ax.hist(
                            pltdf,
                            bins = np.arange(xmin, xmax+(xmax-xmin)/xnbin, (xmax-xmin)/xnbin),
                            weights = weights,
                            density = False,
                            histtype="step",
                            edgecolor = flx_col,
                            #alpha = 0.2,
                            linestyle = "-",
                            linewidth = 2,
                            label = label
                        )
        universes.append(univ)
        
    #-----------------------------------------------------------------#
    #Central Value
    ax.step(bins, which_dict['cv_plot']
         , color = cv_col
         , label =  "Central Value"
    )
    #-----------------------------------------------------------------#
    ax.legend(loc = 'upper left',fontsize = 14)

    plot_tick(ax, 16)
    plot_title(ax, "", 'Opt0 Time Corrected Z % 18.936 [ns]',  "Slices (No Scaling)", 16)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    #-----------------------------------------------------------------#
    fig.tight_layout()

    plt.savefig(savePath+str(name + '_'+which_type + "_weight_universe.png"), dpi=200)
    plt.show()
    
    return universes

#------------------------------------------------------------------------------------------------------------------# 
def loopy_loop_multisim_universe(weight_list, weight_name, len_univ, df, which_dict, which_type, savePath):
    
    cov_array = []
    
    for (name, good_name) in zip(weight_list, weight_name):
        print(name, good_name)
    
        #make new columns by exploding it
        df = make_df_multisim(df, name)
    
        #like what the function says
        var_universe = plot_and_save_universe(df, which_dict
                                              , name, good_name, len_univ
                                              , xmin, xmax, xnbin, which_type
                                              , savePath
                                              )
        #drop columns
        df = df.loc[:,~df.columns.str.startswith(name+'_')]
        
        #make universe into array
        var_universe = np.array(var_universe)
        var_universe = np.sort(var_universe, axis = 0)
    
        #compute cov matrix
        var_cov, _ = get_cov_corr_matrix(which_dict['cv'], var_universe)
    
        #save covariance per variables
        cov_array.append(var_cov)

    return cov_array
#------------------------------------------------------------------------------------------------------------------#   
def check_unisim(row):
    return row[0]

#------------------------------------------------------------------------------------------------------------------# 
def plot_and_save_universe_unisim(df, which_dict, name, good_name, xmin, xmax, xnbin , which_type, savePath):

    ymin, ymax = 0,0
    #-----------------------------------------------------------------#
    if which_type == 'rockbox':
        ymin = rockbox_ymin
        ymax = rockbox_ymax
    elif which_type == 'ncpi0':
        ymin = ncpi0_ymin
        ymax = ncpi0_ymax
    #-----------------------------------------------------------------#
    fig, ax = plt.subplots(1,1, figsize = (6,4))

    #-----------------------------------------------------------------#
    #Universe
    pltdf = df['mod_t']
        
    weights = df[name]
    
    print(len(weights))
    print(len(pltdf))
    
    label = good_name
    
    universe, _, _ = ax.hist(
                            pltdf,
                            bins = np.arange(xmin, xmax+(xmax-xmin)/xnbin, (xmax-xmin)/xnbin),
                            weights = weights,
                            density = False,
                            histtype="step",
                            edgecolor = flx_col,
                            #alpha = 0.2,
                            linestyle = "-",
                            linewidth = 2,
                            label = label
                        )
        
    #-----------------------------------------------------------------#
    #Central Value
    ax.step(bins, which_dict['cv_plot']
         , color = cv_col
         , label =  "Central Value"
    )
    #-----------------------------------------------------------------#
    ax.legend(loc = 'upper left',fontsize = 14)

    plot_tick(ax, 16)
    plot_title(ax, "", 'Opt0 Time Corrected Z % 18.936 [ns]',  "Slices (No Scaling)", 16)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    #-----------------------------------------------------------------#
    fig.tight_layout()

    plt.savefig(savePath+str(name + '_'+ which_type + "_weight_universe.png"), dpi=200)
    plt.show()
    
    return universe
#------------------------------------------------------------------------------------------------------------------# 
def loopy_loop_unisim(df, which_dict, which_type, savePath):
    unisim_cov_array = []
    for name in unisim_list:
        print(name)
    
        universe = plot_and_save_universe_unisim(df, which_dict, name, name, xmin, xmax, xnbin , which_type, savePath)
    
        diff = which_dict['cv'] - universe 
        unisim_cov = np.outer(diff, diff)
        unisim_cov_array.append(unisim_cov)
    
        unisim_err = np.sqrt(np.diag(unisim_cov))
        print(unisim_err)
    
    return unisim_cov_array

#------------------------------------------------------------------------------------------------------------------# 
def throw_universe_tspline(spline, random):
    
    if len(spline) == 0:
        spline = [1, 1, 1, 1, 1, 1]
        
    xs = [-1, 1, -2, 2, -3, 3]
    ys = spline
    
    tsp = ROOT.TSpline3("myspline", array.array('d', xs), array.array('d', ys), len(xs))
    weight = tsp.Eval(random)

    return weight

#------------------------------------------------------------------------------------------------------------------#     
def loopy_loop_multisigma(weight_list, weight_name, random_arr, df, which_dict, which_type, savePath):
    
    cov_array = []
    
    for (name, good_name) in zip(weight_list, weight_name):
        print(name, good_name)
        
        len_univ = np.arange(0, 500)
    
        for idx, r in zip(len_univ, random_arr):
            df[name+'_{}'.format(idx)] = df[name].apply(lambda row: throw_universe_tspline(row, r))
        
        #like what the function says
        var_universe = plot_and_save_universe(df, which_dict
                                              , name, good_name, len(random_arr)
                                              , xmin, xmax, xnbin, which_type
                                              , savePath  
                                              )
        #drop columns
        df = df.loc[:,~df.columns.str.startswith(name+'_')]
        
        #make universe into array
        var_universe = np.array(var_universe)
        var_universe = np.sort(var_universe, axis = 0)
    
        #compute cov matrix
        var_cov, _ = get_cov_corr_matrix(which_dict['cv'], var_universe)
    
        #save covariance per variables
        cov_array.append(var_cov)
        #break

    return cov_array

#------------------------------------------------------------------------------------------------------------------#  
def combine_error(which_dict, error_list):
    
    which_dict['combined_cov'] = new_19_by_19_cov()
    
    for err in error_list:
        which_dict[err +'_cov_frac'] = which_dict[err +'_cov'] / np.outer(which_dict['cv'], which_dict['cv'])
        which_dict[err + '_frac_err'] = np.sqrt(np.diag(which_dict[err +'_cov_frac']))

        which_dict['combined_cov'] = which_dict['combined_cov'] + which_dict[err+'_cov']
        
    which_dict['combined_err'] = np.sqrt(np.diag( which_dict['combined_cov']))
    
    which_dict['combined_cov_frac'] = which_dict['combined_cov'] / np.outer(which_dict['cv'], which_dict['cv'])
    which_dict['combined_frac_err'] = np.sqrt(np.diag(which_dict['combined_cov_frac']))
    
#------------------------------------------------------------------------------------------------------------------#  
def scale_cov_matrix(which_dict, scale_factor, error_list):
    which_dict['cv_scale'] =  which_dict['cv'] * scale_factor
    which_dict['cv_plot_scale'] = which_dict['cv_plot'] * scale_factor
    
    which_dict['combined_cov_scale'] = new_19_by_19_cov()
    
    for err in error_list:
        which_dict[err+'_cov_scale'] = np.outer(which_dict['cv_scale'], which_dict['cv_scale']) * which_dict[err+'_cov_frac']
        which_dict[err+'_err_scale'] = np.sqrt(np.diag(which_dict[err+'_cov_scale']))
        
        which_dict[err+'_cov_frac_scale'] = which_dict[err+'_cov_scale'] / np.outer(which_dict['cv_scale'], which_dict['cv_scale'])
        which_dict[err+'_frac_err_scale'] = np.sqrt(np.diag(which_dict[err+'_cov_frac_scale']))
        
        which_dict['combined_cov_scale'] = which_dict['combined_cov_scale'] + which_dict[err+'_cov_scale']
        
    which_dict['combined_err_scale'] = np.sqrt(np.diag( which_dict['combined_cov_scale']))
    
    which_dict['combined_cov_frac_scale'] = which_dict['combined_cov_scale'] / np.outer(which_dict['cv_scale'], which_dict['cv_scale'])
    which_dict['combined_frac_err_scale'] = np.sqrt(np.diag(which_dict['combined_cov_frac_scale']))
    
#------------------------------------------------------------------------------------------------------------------#  
def plot_combine_err(which_dict, which_type, label, error_list, ifScale = False, scaleYmax = 1, suffix = ''):

    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 2]}, figsize = (6, 6), sharex = True)

    #=========================================================================#
    ymin = 0
    if which_type == 'hnl':
        ymin = hnl_ymin
        ymax = hnl_ymax * scaleYmax
    elif which_type == 'rockbox':
        ymin = rockbox_ymin
        ymax = rockbox_ymax * scaleYmax
    elif which_type == 'ncpi0':
        ymin = ncpi0_ymin
        ymax = ncpi0_ymax * scaleYmax
    elif which_type == 'nu':
        ymin = nu_ymin
        ymax = nu_ymax * scaleYmax
    elif which_type == 'cos':
        ymin = cos_ymin
        ymax = cos_ymax * scaleYmax
    
    #=========================================================================#
    #central value
    ax1.step(bins, which_dict['cv_plot' + suffix]
         , color = cv_col
         , label =  label
        )
    #-----------------------------------------------------------------#
    #universe 1 sigma
    
    prefix = 'combined_err'
    col = total_col
    label = 'Total'
    if which_type == 'cos':
        prefix = 'stat_err'
        col = stats_col
        label = "Satistics"

    bottom = which_dict['cv' + suffix] - which_dict[prefix + suffix]
    peak = which_dict['cv' + suffix] + which_dict[prefix + suffix]
    bottom = np.nan_to_num(bottom)
    peak = np.nan_to_num(peak)
    height = peak - bottom

    plt. rcParams["hatch.color"] = col
    ax1.bar(
        x = bins[:-1]
        , height= height
        , width = np.diff(bins)
        , bottom = bottom
        , align = 'edge'
        , hatch='///'
        , fc=col
        , alpha = 0.2
        , label = label
                        )

    #-----------------------------------------------------------------#
    plot_tick(ax1, 16)
    
    ytitle =  r"Slices (No Scaling)"
    if ifScale == True:
         ytitle = r"Slices (1$\times10^{21}$ POT)"
            
    plot_title(ax1, "", '',  ytitle, 16)

    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)

    ax1.legend(loc='upper left', fontsize=14, fancybox=False, ncol = 1)
    
    #=========================================================================#
    
    #-----------------------------------------------------------------#
    for err in error_list:
        
        if err == 'stat':
            error_label = 'Statistics'
            col = stats_col
        elif err == 'flx':
            error_label = 'Flux'
            col = flx_col    
        elif err == 'xsec':
            error_label = 'Cross Section'
            col = xsec_col
        elif err == 'g4':
            error_label = 'Re-Interaction'
            col = g4_col
            
        ax2.step(bins, np.insert(which_dict[err + '_frac_err' + suffix], 0,0) 
         , color = col , linestyle = '-'
         , label = error_label, lw = 2
        )
    #-----------------------------------------------------------------#
    #COMBINED
    if which_type != 'cos':

        ax2.step(bins, np.insert(which_dict['combined_frac_err' + suffix], 0,0) 
         , color = total_col , linestyle = '--'
         , label = 'Total', lw = 2
        )
    #-----------------------------------------------------------------#
    if which_type == 'hnl':
        ymax = hnl_ymax2
    elif which_type == 'rockbox':
        ymax = rockbox_ymax2
    elif which_type == 'ncpi0':
        ymax = ncpi0_ymax2
    elif which_type == 'nu':
        ymax = nu_ymax2
    elif which_type == 'cos':
        ymax = cos_ymax2
    
    plot_tick(ax2, 16)
    plot_title(ax2,"", 'Opt0 Time Corrected Z % 18.936 [ns]', "Fractional Error", 16)


    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)
    ax2.legend(loc='upper center', fontsize=12, fancybox=False, ncol = 2)
    #=========================================================================#
    fig.tight_layout()
