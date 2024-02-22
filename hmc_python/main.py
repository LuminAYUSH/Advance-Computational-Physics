import numpy as np
import scipy
import sklearn as sk
import matplotlib.pyplot as plt
import sys
import os
import math

from util import *

PI = math.pi
COLDLAT = 0
HOTLAT = 1
PLUS = 1
MINUS = -1
EVEN = 0x02
ODD = 0x01
EVENANDODD = 0x03
LSIZE = 16   # Change LSIZE to 16/32 accordingly
DELTAMAX = 50.0

# acl parameters for switch = 0

T_cut = 5
D_cut = 5
MAXT_cut = 25
NOT_cut = 5

#------------------------------------------------
#Declared external variables in the header file
#------------------------------------------------

nx = None
nt = None
volume = None
nf = None
mdstep = None
cgiter1 = None
cgiter2 = None
# long *iseed # Check this later !!!!
g = None
step = None
residue1 = None
residue2 = None
mid = None
no_even_sites = None
no_odd_sites = None
no_garbage = None
bin_length = None
no_bin = None
meas_loop = None
meas_length = None
prop_length = None
no_meas = None
seg_length = None
no_a_seg = None
no_prop_seg = None
hmc_it = None
counter = None
sw_flag = None #flag to switch action between garbage and
# #                  autocorln. loops and measurement loops.
# #                  sw_flag = 0 => garbage & autocorln. loop;
# #                  sw_flag = 1 => measurement loop; */

lattice = None
store = None
conf = None
con = None
garbage = None
ac_store = None
ac_prop = None
bin_av = None
psi = None
psi_acl = None
G_prop = None
G_store = None
G_temp = None
prop = None
tprop = None
T_int = None
T_int_prop = None
gen_pt = None

#------------------------------------------------
# File IO Defined in the code
#------------------------------------------------
pin = open("sigma1.in","r")
ptout = open("sigma1.out", "a")
ptacl = open("sigma1.acl", "a");
ptlat = open("sigma1.lat", "a");
ptprop = open("sigma1.prop", "w");
ptpropacl = open("sigma1.propacl", "w");

prompt = setup_gn()
readin(prompt)

# randomize() # Not defined here but present in the original code

# filelat(pin) # Getting converted will be merged later

av_sigma = 0.0  # grand average of <sigma>
av_psi = 0.0  # grand average of <psi_bar-psi>
t_ex_sigma = 0.0  # grand sum over all <sigma>s
seg_av_sigma = 0.0  # segment-average of <sigma> for tau_int
seg_av_prop = 0.0  # segment-average of <propagator[m]>

no_hmc = 0  # # of hmc steps calculated
no_acc = 0  # # of accepted configurations
counter = 0  # counter used in hmc.c
meas = 0  # meas is a switch for measurements
no_auto = 0  # no_auto is the index for autocoreln. measurement
no_prop = 0  # no_prop is the index for Prop-Acl calculation

a = 0  # a index used in tau_int measurements
ac = 0  # ac configuration index to ac_prop[ac]
acc = 0  # acc index used in G_temp[][acc] i.e. # of data point and also used as index in ac_store
bin = 0  # bin is index to bin_average
g = 0  # g configuration index to garbage i.e. garbage[g]
k = 0  # k configuration index to store i.e. store[k]
j = 0  # j configuration index to store[j], during measurements

if sw_flag == 0:
    #  GARBAGE LOOPS AND AUTOCORELATION MEASUREMENTS LOOPS
    for n in range(hmc_it):
        no_acc += 1
        no_hmc += hmc()

        if no_acc<=no_garbage:
            garbage[g] = average_sigma()
            g += 1
            bin +=1
            if bin%bin_length == 0:
                for m in range((g-bin),g):
                    bin_av[k] +=garbage[m]/bin_length
                k += 1
                bin = 0

        if no_acc>no_garbage:
            ac_store[acc] = average_sigma()
            acc += 1
            no_auto += 1

            if no_auto%seg_length == 0:
                lbd = acc - no_auto
                for m in range(lbd,acc):
                    seg_av_sigma += ac_store[m]/seg_length
                autocorel(seg_av_sigma,lbd,a)
                a += 1
                seg_av_sigma = 0
                no_auto = 0

        for i in range(0,volume):
            con[i] = conf[i]
        for i in range(0,volume):
            conf[i] = lattice[i].sigma

    acc_rate = no_acc/no_hmc

    for n in range(0, MAXT_cut):
        ptacl.write(f"{n+1}")
        for u in range(NOT_cut):
            d_t_int = 0
            av_t_int = 0

            for a in range(no_a_seg):
                av_t_int += T_int[u][n][a]/no_a_seg
            for a in range(no_a_seg):
                d_t_int += (T_int[u][n][a] - av_t_int)**2

            d_t_int = math.sqrt(d_t_int/(no_a_seg-1))

            ptacl.write(f"\t{av_t_int}\t{d_t_int}")
        ptacl.write("\n")

    print(f"\n\n no_acc_traj={no_acc} \t no_hmc={no_hmc} \t acc_rate={acc_rate}")

    for k in range(0,no_bin):
        ptout.write(f"{k+1}\t{bin_av[k]}\n")
    for i in range(0,volume):
        ptlat.write(f"{lattice[i].sigma}\n")



if sw_flag == 1:
    for n in range(meas_loop):
        no_acc += 1
        no_hmc += hmc()

        meas = meas + 1
        no_prop = no_prop + 1

        if meas%meas_length == 0:
            store[j] += average_sigma()
            t_ex_sigma += store[j]
            j += 1

        if meas%prop_length == 0:
            pass

    av_sigma = t_ex_sigma / no_meas
    sqdev = 0
    for k in range(no_meas):
        sqdev += (store[k]-av_sigma)**2
    d_av_sigma = sqrt(sqdev/(no_meas-1))

    acc_rate = no_acc/no_hmc

    print(f"\n\n\ no_acc_traj = {no_acc} \t no_hmc = {no_hmc} \t acc_rate = {acc_rate} \n\n")
    print(f"av_sigma={av_sigma} \t av_psi{av_psi}\n\n")
    print(f"d_av_sigma={d_av_sigma} \t d_av_psi{d_av_psi}\n\n")

    for i in range(0,volume):
        ptlat.write(f"{lattice[i].sigma}\n")

else:
    print("KILL YOURSELF, PLEASE !!!")
