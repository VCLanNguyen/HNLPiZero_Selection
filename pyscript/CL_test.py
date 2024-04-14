import sys
import os

os.nice(20)

import warnings
warnings.filterwarnings("ignore")

import argparse                                   
import numpy as np
import math
import pyhf
pyhf.set_backend("numpy")

# Local helper script
hnlDIR = os.environ['_']
sys.path.append(hnlDIR + '/pyscript/')

#------------------------------------------------------------------------------------------------------------------#
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()

    #print("Found nearest value = " + str(array[idx]))
    #print("Nearest value idx = " + str(idx))
    return idx

#------------------------------------------------------------------------------------------------------------------#
def find_crit_hypo_test(signal, bkg, sigma_bkg, obs):

    model = pyhf.simplemodels.uncorrelated_background(
            signal=[signal], 
            bkg=[bkg], 
            bkg_uncertainty=[sigma_bkg]
            )

    observations = [obs]

    data = pyhf.tensorlib.astensor(observations + model.config.auxdata)

    mu_test = 1.0

    CLs_obs, CLs_exp, CLs_exp_band = pyhf.infer.hypotest(
                            mu_test, data, model,
                            return_expected=True,
                            return_expected_set=True, 
                            test_stat="qtilde"
                            )

    #print("Confidence Level Observed")
    #print(CLs_obs)
    #print("Confidence Level Expected")
    #print(CLs_exp)
    #print("Confidence Level Expected Band")
    #print(CLs_exp_band)

    return CLs_obs, CLs_exp, CLs_exp_band
#------------------------------------------------------------------------------------------------------------------#
def find_90CLs(Umu, signal, bkg, sigma_bkg, obs):

    Umu_arr = np.logspace(-8, -6, 100, endpoint = True)

    signal_arr = []

    for u in Umu_arr:
        s = signal * (u/ Umu)**2
        signal_arr.append(s)


    cl_obs_arr = []
    cl_exp_arr = []
    cl_band_arr = []

    for s in signal_arr:
        cl_obs, cl_exp, cl_band = find_crit_hypo_test(s, bkg, sigma_bkg, bkg)

        cl_obs_arr.append(cl_obs)
        cl_exp_arr.append(cl_exp)
        cl_band_arr.append(cl_band)

    Umu_arr = np.array(Umu_arr)
    signal_arr = np.array(signal_arr)
    cl_obs_arr = np.array(cl_obs_arr)
    cl_exp_arr = np.array(cl_exp_arr)
    cl_band_arr = np.array(cl_band_arr)
    
#    for u, s, cl_obs, cl_exp in zip(Umu_arr, signal_arr, cl_obs_arr, cl_exp_arr):
#        print("Umu = {0:.9f}, nsignal = {1:.3g}, CLs_obs = {2:.3g}, CLs_exp = {3:.3g}".format(u, s, cl_obs, cl_exp))

    #90% CLs --> CLs = 0.01
    idx = find_nearest_idx(cl_obs_arr, value = 0.01)
    Umu_crit = Umu_arr[idx]
    signal_crit = signal_arr[idx]
    cl_obs_crit = cl_obs_arr[idx]
    cl_exp_crit = cl_exp_arr[idx]
    cl_band_crit = cl_band_arr[idx]

    print("Critical coupling = " + str(Umu_crit))
    print("Critical nSignals = " + str(signal_crit))
    print("Critical CLs obs = " + str(cl_obs_crit))
    print("Critical CLs exp = " + str(cl_exp_crit))
    print("Critical CLs band = " + str(cl_band_crit))
    print("\n")

    for cl in cl_band_crit:
        idx = find_nearest_idx(cl_exp_arr, value = cl)
        print("CL = " +str(cl))
        print("Umu = " +str(Umu_arr[idx]))
        print("nSignals = " +str(signal_arr[idx]))
        print("\n")

#------------------------------------------------------------------------------------------------------------------#
def main(args):                                            

    Umu = 1e-7

    signal = 45
    bkg = 1
    sigma_bkg = np.sqrt(bkg)

    print("--------------------------------------------")
    print("Mixing Angle = " + str(Umu))
    print("Signal Events = " + str(signal))
    print("Bkg Events = " + str(bkg))
    print("Sig/Bkg = " + str(signal/bkg*100))
    print("Bkg Error = " + str(sigma_bkg))
    print("--------------------------------------------")

    find_90CLs(Umu, signal, bkg, sigma_bkg, bkg)



if __name__ == '__main__':                                 
                                                           
    parser = argparse.ArgumentParser()                     
    args = parser.parse_args()                             
    main(args)                                             
