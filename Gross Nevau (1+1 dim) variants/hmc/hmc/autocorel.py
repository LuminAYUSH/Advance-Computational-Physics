#/* ********** autocorel.c ********** */
#/* *** Alpha Machine, version 2 ***  */
#/* This routine calculates the unnormalized autocorrelated functions and
#   normalized autocorrelated function along with integrated autocorrela-
#   tion time. It reads some T_cut and calculates var(int. ac time).  */
   
from lattice_gn import T_cut, MAX_Tcut, seg_length, ac_store, NOT_cut, rho, T_int, D_cut
import numpy as np

def autocorel(sigma_av, lb, a_index):
    
    global T_cut, MAX_Tcut, seg_length, ac_store, NOT_cut, rho, T_int, D_cut
    
    c0 = 0.0; N_t = 0; tcut = T_cut
    
    #Calculating the unnormalized autocorrelation functions ct and c0
    #and normalized autocorrelation function rho.
    
    for j in range(lb, lb+seg_length):
        c0 += ((ac_store[j] - sigma_av) * (ac_store[j] - sigma_av))/seg_length
    
    for u in range(0, NOT_cut):
        for t in range(1, tcut+1):
            N_t = seg_length - t
            ct = 0.0
            for j in range(lb, lb+N_t):
                ct += ((ac_store[j] - sigma_av) * (ac_store[j+t] - sigma_av))/N_t
                rho[t-1] += ct/c0
        
        for t in range(0, tcut): T_int[u][t][a_index] = 0.5
        
        #Calculation of tau_int, the integrated autocorrelation time
        for t in range(0, tcut):
            for k in range(0, t+1): T_int[u][t][a_index] += rho[k]
        
        tcut += D_cut
        
def average_sigma(sites):
    t_sigma = sum(site.sigma for site in sites)
    av_sigma = t_sigma / len(sites)
    return av_sigma