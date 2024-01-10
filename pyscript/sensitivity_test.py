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
def Asimov_Sensitivity_Sigma(s, b, sigb):

    part1 =  ((s+b) * (b + sigb**2)) / (b**2 + (s+b) * sigb**2)
    part2 = 1 + (sigb**2 / ( b * (b + sigb**2)))
    lnpart1 = math.log(part1)
    lnpart2 = math.log(part2)

    sigma = np.sqrt(2 * ( (s+b)*lnpart1 - (b**2/sigb**2)*lnpart2 ))

    return sigma
#------------------------------------------------------------------------------------------------------------------#
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()

    #print("Found nearest value = " + str(array[idx]))
    #print("Nearest value idx = " + str(idx))
    return idx
#------------------------------------------------------------------------------------------------------------------#
def find_crit_Umu_Asimov(Umu, signal, bkg, sigma_bkg):

    print("Asimov")

    Umu_arr = np.logspace(-9, -7, 1000000, endpoint = True)

    signal_arr = []

    for u in Umu_arr:
        s = signal * (u/ Umu)**2
        signal_arr.append(s)


    ss_arr = []

    for s in signal_arr:
        ss = Asimov_Sensitivity_Sigma(s, bkg, sigma_bkg)
        ss_arr.append(ss)
    
    Umu_arr = np.array(Umu_arr)
    signal_arr = np.array(signal_arr)
    ss_arr = np.array(ss_arr)

    #drop nan values from the asimov equation
    mask = ~np.isnan(ss_arr)
    Umu_arr = Umu_arr[mask]
    signal_arr = signal_arr[mask]
    ss_arr = ss_arr[mask]

    #90% CL = 1.645 sigma from the mean
    idx = find_nearest_idx(ss_arr, value = 1.645)
    Umu_crit = Umu_arr[idx]
    ss_crit = ss_arr[idx]

    print("--------------------------------------------")
    print("Crticical sensitivity sigma = " + str(ss_crit))
    print("Critical Umu = " + str(Umu_crit))
    print("--------------------------------------------")

    return Umu_crit

#------------------------------------------------------------------------------------------------------------------#
def find_crit_Umu_Gaus(Umu, signal, bkg):
    
    print("Folded Gaussian")
    
    Umu_arr = np.logspace(-9, -7, 1000000, endpoint = True)

    signal_arr = []

    for u in Umu_arr:
        s = signal * (u/ Umu)**2
        signal_arr.append(s)

    crit_val = 1.282 * np.sqrt(bkg)

    idx = find_nearest_idx(signal_arr, value = crit_val)
    Umu_crit = Umu_arr[idx]
    signal_crit = signal_arr[idx]

    print("--------------------------------------------")
    print("Crticical signal events = " + str(signal_crit))
    print("Critical Umu = " + str(Umu_crit))
    print("--------------------------------------------")

    return Umu_crit
    
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

    CLs_obs, CLs_exp_band = pyhf.infer.hypotest(
    mu_test, data, model, return_expected_set=True, test_stat="qtilde"
    )

    print(CLs_obs)
    print(CLs_exp_band)
#------------------------------------------------------------------------------------------------------------------#
def main(args):                                            

    print("5 ns SELECTION")
    Umu = 1e-7
    signal = 75
    bkg = 184
    sigma_bkg = np.sqrt(bkg)

    print("--------------------------------------------")
    print("Mixing Angle = " + str(Umu))
    print("Signal Events = " + str(signal))
    print("Bkg Events = " + str(bkg))
    print("Sig/Bkg = " + str(signal/bkg*100))
    print("Bkg Error = " + str(sigma_bkg))
    print("--------------------------------------------")


    find_crit_hypo_test(signal, bkg, sigma_bkg, bkg)
    find_crit_hypo_test(signal, bkg, sigma_bkg, signal+bkg)

    print("3 ns SELECTION")
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


    find_crit_hypo_test(signal, bkg, sigma_bkg, bkg)
    find_crit_hypo_test(signal, bkg, sigma_bkg, signal+bkg)


if __name__ == '__main__':                                 
                                                           
    parser = argparse.ArgumentParser()                     
    args = parser.parse_args()                             
    main(args)                                             
