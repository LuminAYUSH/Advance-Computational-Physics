### Util for the HMC code with numpy vectorisation :>
import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys
import os
import math
import random
from tqdm.auto import tqdm
from .lparam import *
import numpy as np



class LATTICE:
    def __init__(self,volume_):
        self.x = np.array(list(None for i in range(volume_)))
        self.t = np.array(list(None for i in range(volume_)))
        self.sign = np.array(list(None for i in range(volume_)))
        self.parity = np.array(list(None for i in range(volume_)))
        self.sigma = np.array(list(None for i in range(volume_)))
        self.phi = np.array(list(None for i in range(volume_)))
        self.mom = np.array(list(None for i in range(volume_)))
        self.chi = np.array(list(list(None for i in range(4)) for j in range(volume_)))
        self.eta = np.array(list(list(None for i in range(4)) for j in range(volume_)))
        self.p = np.array(list(None for i in range(volume_)))
        self.r = np.array(list(None for i in range(volume_)))
        self.mp = np.array(list(None for i in range(volume_)))
        self.mmp = np.array(list(None for i in range(volume_)))

    def __getitem__(self,key):
        return self.__dict__[key]
    def __setitem__(self,key,value):
        self.__dict__[key] = value


def setup_gn():
    prompt = initial_set()
    layout()
    make_lattice()
    make_nn_gather()
    return prompt


def initial_set():
    global nx
    global nt
    global sw_flag
    global no_garbage
    global bin_length
    global hmc_it
    global meas_length
    global prop_length
    global seg_length
    global volume
    global mid
    global meas_loop

    prompt = None
    status = None

    print("GN model with HMC algorithm\n")
    print("Alpha machine, Version 1\n")
    print("type 0 for no prompts, 1 for prompts or 2 for list of prompts\n")

    try:
        prompt = getprompt() #Using the direct value inplace of the pointer in the original code base
    except:
        print("error in input: initial prompt")
        return -1

#     nx = get_i(prompt, "nx")
    nx = 16
#     nt = get_i(prompt, "nt")
    nt = 16

    if nx%2 !=0 and nt%2 !=0:
        print("nx, nt must be even!! \n")
        sys.exit()

    print(f"Lattice dimensions = {nx} {nt}\n")

#     Switch flag
#     sw_flag = get_i(prompt,"switch_flag")
    sw_flag = 0
    print(f"Switch_Flag = {sw_flag}\n")

#     Number of measurements
#     no_garbage = get_i(prompt,"no_of_garbage_loops")
    no_garbage = 100
    print(f"No of garbage loops = {no_garbage}")

    # the length of each bin
#     bin_length = get_i(prompt,"bin_length")
    bin_length = 5
    print(f"bin_length = {bin_length}")

    # Number of HMC iterations, only accepted ones count
#     hmc_it = get_i(prompt,"no_of_hmc_iterations")
    hmc_it = 1100
    print(f"# of hmc iterations = {hmc_it}")

    # length after which fermionic observables are measured
#     meas_length = get_i(prompt,"meas_length")
    meas_length = 5
    print(f"meas_length = {meas_length}")

    # length after which the propagator autocorrelations are calculated
#     prop_length = get_i(prompt,"prop_length")
    prop_length = 5
    print(f"prop_length = {prop_length}")

    # segment length after which autocorrelations are measured
#     seg_length = get_i(prompt,"seg_length")
    seg_length = 100
    print(f"seg_length = {seg_length}")

    volume = nx*nt
    mid = nt/2
    meas_loop = hmc_it-no_garbage
    if sw_flag == 1:
        seg_length = prop_length
    return prompt


def make_lattice():
    global no_bin
    global no_meas
    global no_a_seg
    global no_prop_seg
    global no_garbage
    global meas_length
    global meas_loop
    global seg_length
    global prop_length
    global volume
    global nt
    global MAXT_cut
    global NOT_cut

    global lattice
    global store
    global conf
    global con
    global garbage
    global ac_store
    global ac_prop
    global bin_av
    global psi
    global psi_acl
    global G_prop
    global G_store
    global G_temp
    global prop
    global tprop
    global T_int
    global T_int_prop
    global neighbor
    global gen_pt

    i = None
    j = None
    k = None

    x = None
    t = None

    no_bin = int(no_garbage/bin_length)
    no_meas = int(meas_loop/meas_length)
    no_a_seg = int(meas_loop/seg_length)
    no_prop_seg = int(meas_loop/prop_length)

    #All allocations managed by python backend
    lattice = LATTICE(volume)
    
    store = np.array(list(None for i in range(no_meas)))
    conf = np.array(list(None for i in range(volume)))
    con = np.array(list(None for i in range(volume)))
    garbage = np.array(list(None for i in range(no_garbage)))
    ac_store = np.array(list(None for i in range(meas_loop)))
    ac_prop = np.array(list(None for i in range(prop_length)))
    bin_av = np.array(list(None for i in range(no_bin)))
    psi = np.array(list(None for i in range(no_meas)))
    psi_acl = np.array(list(None for i in range(meas_loop)))
    G_prop = np.array(list(None for i in range(nt)))
    G_store = np.array(list(list(None for i in range(no_meas)) for j in range(nt))) # This should be nx see the defination using LSIZE
    G_temp = np.array(list(list(None for i in range(meas_loop)) for j in range(nt)))  # This should be nx see the defination using LSIZE
    prop = np.array(list(None for i in range(nt)))
    tprop = np.array(list(None for i in range(nt)))

    T_int = np.array(list(list(list(None for i in range(no_a_seg)) for j in range(MAXT_cut)) for k in range(NOT_cut)))
    T_int_prop = np.array(list(list(list(None for i in range(no_prop_seg)) for j in range(MAXT_cut)) for k in range(NOT_cut)))
    neighbor = np.array(list(list(None for i in range(volume)) for j in range(4)))
    gen_pt = np.array(list(None for i in range(volume))) # Not vectorised for sake of simplicity

#################################################
#Vectorise Later
#################################################

    for t_ in range(nt):
        for x_ in range(nx):
            i = site_index(x_,t_)     # Function not defined yet !!!!
            lattice.x[i] = x_
            lattice.t[i] = t_
            if t_%2 ==0:
                lattice.sign[i] = 1
            else:
                lattice.sign[i] = -1
                if (x_+t_)%2 == 0:
                    lattice.parity[i] = EVEN
                else:
                    lattice.parity[i] = ODD
                    
#################################################
#################################################

    ac_store[:] = 0
    bin_av[:] = 0
    G_prop[:] = 0
    T_int[:,:,:] = 0


def get_f(prompt,variable_name_string):

    if prompt == 1:
        x = float(input(f"enter {variable_name_string}"))
        return(x)

    else:
        print("Change Prompt type currently not supported")
        system.exit()


def get_i(prompt, variable_name_string):

    if prompt == 1:
        x = int(input(f"enter {variable_name_string}"))
        return(x)

    else:
        print("Change Prompt type currently not supported")
        system.exit()


def getprompt():
    prompt = int(input("Enter the prompt type (only 1 is supported for now):"))
    return prompt


def readin(prompt):

    global g
    global nf
    global mdstep
    global step
    global cgiter1, cgiter2
    global residue1, residue2

    status = None
    x = None

#     g = get_f(prompt,"g")
    g = 0.46
#     nf = get_i(prompt,"nf")
    nf = 4

#     mdstep = get_i(prompt,"no_of_md_steps")
    mdstep = 25
#     step = get_f(prompt,"step_size")
    step = 0.04
    print(f"no of md steps = {mdstep}, step size = {step}")

#     cgiter1 = get_i(prompt,"max_cg_iterations_for_hamil")
#     cgiter2 = get_i(prompt,"max_cg_iterations_for_piup")
    cgiter1 = 5000
    cgiter2 = 2000
    print(f"maximum no. of conj. grad. iterations = {cgiter1},{cgiter2}")

#     residue1=get_f(prompt,"residue_for_cg_hamil")
#     residue2=get_f(prompt,"residue_for_cg_piup")
    residue1 = 1e-7
    residue2 = 1e-5
    print(f"residues for conjugate grad = {residue1},{residue2}")


def autocorel(sigma_av, lb, a_index):

    global T_cut, MAX_Tcut, seg_length, ac_store, NOT_cut, T_int, D_cut
    rho = np.zeros(MAX_Tcut)

    c0 = 0.0; N_t = 0; tcut = T_cut


    #Calculating the unnormalized autocorrelation functions ct and c0
    #and normalized autocorrelation function rho.        
    c0 = np.sum((ac_store[lb:(lb+seg_length)] - sigma_av)**2)/seg_length

    
    for u in range(0, NOT_cut):
        for t in range(1, tcut+1):
            N_t = seg_length - t
            ct = np.sum(((ac_store[lb:(lb+N_t)] - sigma_av)**2)/N_t)
            rho[t-1] += ct/c0

        T_int[u][:tcut][a_index] = 0.5

        #Calculation of tau_int, the integrated autocorrelation time
        for t in range(0, tcut):
                T_int[u][t][a_index] = np.sum(rho[:(t+1)])

        tcut += D_cut


def average_sigma():
    global volume, lattice
    i = None
    av_sigma = None
    t_sigma = None
    t_sigma = 0
    
    t_sigma = np.sum(lattice.sigma[:])
    av_sigma = t_sigma/volume
    return av_sigma



def cg_md(src,dest,cgiter,residue,cgflag,flavor):
    N_iter = None
    i = None
    size_src = None
    cp = None
    d = None
    dsize_r = None
    dsize_src = None
    size_r = None
    a = None
    b = None
    c = None

    global volume
    global PLUS
    global EVENANDODD
    global lattice

    dsize_src = 0
    
    dsize_src = np.sum((lattice[src][:,flavor])**2)
    size_src = math.sqrt(dsize_src)

    if cgflag == 0:
        
        lattice[dest][:,flavor] = 0
        lattice.r[:] = lattice[src][:,flavor]
        lattice.p[:] = lattice.r[:]
                       
        dsize_r = 1
        size_r = dsize_r

    if cgflag != 0:
        matp2d(dest,"p",PLUS,EVENANDODD,flavor)
        matd2d("p","mp",MINUS,EVENANDODD)

        dsize_r = 0

        lattice.r[:] = lattice[src][:,flavor] - lattice.mp[:]
        dsize_r = np.sum(lattice.r[:]**2)
        lattice.p[:] = lattice.r[:]

        size_r = math.sqrt(dsize_r)/size_src

    cp = dsize_r

    N_iter = 0

    while N_iter < cgiter and size_r > residue:
        c = cp

        matd2d("p","mp",PLUS,EVENANDODD)

        d=0
        d = np.sum((lattice.mp[:])**2)
        a = c/d

        matd2d("mp","mmp",MINUS,EVENANDODD)

        cp = 0
        lattice[dest][:,flavor] = lattice[dest][:,flavor] + a*lattice.p[:]
        lattice.r[:] = lattice.r[:] - a*lattice.mmp[:]
        cp = np.sum(lattice.r[:]**2)

        b = cp/c
        dsize_r = 0
                       
        lattice.p[:] = lattice.r[:] + b*lattice.p[:]
        dsize_r = np.sum(lattice.r[i]**2)
        size_r = math.sqrt(dsize_r)/size_src

        N_iter = N_iter + 1


    if size_r > residue:
        print("CG_MD Not Converged")
        system.exit(1)


def cg_prop(src , dest , cgiter , residue, cgflag):
    global volume, lattice, PLUS, MINUS, EVENANDODD, lattice
    N_iter = 0
    size_src = 0.0
    cp = 0.0
    dsize_r = 0.0
    dsize_src = 0.0
    size_r = 0.0
    a = 0.0
    b = 0.0
    c = 0.0

    # Normalisation
    dsize_src = np.sum(lattice[src][:]**2)
    size_src = math.sqrt(dsize_src)

    # Initial guess
    if cgflag == 0:
                       
        lattice[dest][:] = 0
        lattice.r[:] = lattice[src][:]
        
        dsize_r = 1.0
        size_r = dsize_r

    if cgflag != 0:
                       
        matd2d(dest,"mp",PLUS,EVENANDODD)
        dsize_r = 0.0
                       
        lattice.r[:] = lattice[src][:] - lattice.mp[:]
        dsize_r = np.sum(lattice.r[:]**2)
        size_r = math.sqrt(dsize_r)/size_src

    matd2d("r","p",MINUS,EVENANDODD)

    cp = 0.0
    
    cp = np.sum(lattice.p[:]**2)
    # Start of CG iteration loop
                       
    while N_iter < cgiter and size_r > residue:
        c = cp

        matd2d("p","mp",PLUS,EVENANDODD)

        d = 0.0
        d = np.sum(lattice.mp[:]**2)
        a = c/d


        lattice[dest][:] = lattice[dest][:] + a*lattice.p[:]
        lattice.r[:] = lattice.r[:] - a*lattice.mp[:]

        matd2d("r","mp",MINUS,EVENANDODD)

        cp = 0.0                       
        cp = np.sum(lattice.mp[:]**2)

        b = cp/c

        dsize_r = 0.0

        lattice.p[:] = lattice.mp[:] + b*lattice.p[:]
        
        dsize_r = np.sum(lattice.r[:]**2)
        size_r = math.sqrt(dsize_r)/size_src
        N_iter += 1

    if size_r > residue:
        print("CG_PROP Not Converged")


def matd2d(src,dest,isign,parity):
    i = None
    n = None

    global XUP,XDN,TUP,TDN, EVENANDODD, volume, lattice

    gather(src, XUP, EVENANDODD, gen_pt)
    lattice[dest][:] = 0.5*lattice.sign[:]*gen_pt[:]

    gather(src, XDN, EVENANDODD, gen_pt)
    lattice[dest][:] = lattice[dest][:] - 0.5*lattice.sign[:] * gen_pt[:]

    gather(src, TUP, EVENANDODD, gen_pt)
    lattice[dest][:] = lattice[dest][:] + 0.5 * gen_pt[:]

    gather(src, TDN, EVENANDODD, gen_pt)
    lattice[dest][:] = lattice[dest][:] - 0.5 * gen_pt[:]

    lattice[dest][:] = isign*(lattice[dest][:]) + (lattice.phi[:])*(lattice[src][:])


def matd2p(src,dest,isign,parity,flavor):
    i = None
    n = None

    global volume, XUP, XDN,gen_pt, TDN, TUP, EVENANDODD, lattice

    gather(src, XUP, EVENANDODD, gen_pt)
    lattice[dest][i,flavor] = lattice.sign[:] * 0.5 * gen_pt[:]

    gather(src, XDN, EVENANDODD, gen_pt)
    lattice[dest][:,flavor] = lattice[dest][:,flavor] - 0.5*lattice.sign[:] * gen_pt[:]
                       
    gather(src, TUP, EVENANDODD, gen_pt)
    lattice[dest][:,flavor] = lattice[dest][:,flavor] + 0.5 * gen_pt[:]

    gather(src, TDN, EVENANDODD, gen_pt)
    lattice[dest][:,flavor] = lattice[dest][:,flavor] - 0.5 * gen_pt[:]

    lattice[dest][:,flavor] = isign*lattice[dest][:,flavor] + lattice.phi[:]*lattice[src][:]


def matp2d(src,dest,isign,parity,flavor):
    i = None
    n = None

    global volume, XUP,gen_pt, XDN, TDN, TUP, EVENANDODD, lattice

    gather(src, XUP, EVENANDODD, gen_pt)
    lattice[dest][:] = 0.5*lattice.sign[:] * gen_pt[:,flavor]

    gather(src, XDN, EVENANDODD, gen_pt)
    lattice[dest][:] = lattice[dest][:] - 0.5 * lattice.sign[:] * gen_pt[:,flavor]

    gather(src, TUP, EVENANDODD, gen_pt)
    lattice[dest][:] =  lattice[dest][:] + 0.5 * gen_pt[:,flavor]

    gather(src, TDN, EVENANDODD, gen_pt)
    lattice[dest][:] = lattice[dest][:] - 0.5 * gen_pt[:,flavor]

    lattice[dest][:] = isign*lattice[dest][:] +lattice.phi[:]*lattice[src][:,flavor]



def matp2p(src,dest,isign,parity,flavor):
    i = None
    n = None

    global XUP,XDN,TUP,TDN,gen_pt,EVENANDODD, volume, lattice
    gather(src, XUP, EVENANDODD, gen_pt)
    lattice[dest][:,flavor] = lattice.sign[:] * 0.5 * gen_pt[:,flavor]

    gather(src, XDN, EVENANDODD, gen_pt)
    lattice[dest][:,flavor] -= lattice.sign[:] * 0.5 * gen_pt[:,flavor]

    gather(src, TUP, EVENANDODD, gen_pt)
    lattice[dest][:,flavor] +=  0.5 * gen_pt[:,flavor]

    gather(src, TDN, EVENANDODD, gen_pt)
    lattice[dest][:,flavor] -=  0.5 * gen_pt[:,flavor]

                       
    lattice[dest][:,flavor] = isign*(lattice[dest][:,flavor]) + (lattice.phi[:])*(lattice[src][:,flavor])
    


def propagator():

    prompt = None
    i = None
    m = None
    n = None
    n0 = None
    t = None
    x = None
    xl = None
    xu = None
    source = None
    lspi = None


    global lattice, volume, prop, tprop, nt, nx, mid

    lpsi = 0.0

    tprop = np.array(list(0.0 for i in range(nt)))

    lattice.r[:] = 0.0

    for t in range(nt):
        source = (t * nx) + mid
        lattice[source].r[:] = 1.0

        cg_prop("r", "mmp", cgiter1, residue1, 0)

        for x in range(nx):
            xl = x * nx
            xu = ((x + 1) * nx) - 1
            n0 = x - t
            if n0 < 0:
                n0 += nt
            prop[n0] = 0.0

            for m in range(xl, xu):
                prop[n0] += lattice.mmp[m]

            tprop[n0] += prop[n0]

        lpsi += lattice.mmp[source]

    tprop = tprop/nt
    lpsi = lpsi / nt

    return lpsi


def make_nn_gather():
    #This will run only once hence not vectorised
    global neighbor, lattice
    i = None
    j = None
    dir = None

    xpt = None
    tpt = None
    mul = None

    for dir in range(0,4):
        for i in range(volume):
            xpt,tpt = neighbor_coords(lattice.x[i],lattice.t[i], dir)
            j = site_index(xpt,tpt)
            neighbor[dir,i] = j
    neighbor = neighbor.astype(np.uint16)

def neighbor_coords(x,t,dir):
    global neighbor
    xp = x
    tp = t

    global nx, nt

    if dir == 0:
        xp = (x+1)%nx
    if dir == 3:
        xp = (x+nx-1)%nx
    if dir == 1:
        tp = (t+1)%nt
    if dir == 2:
        tp = (t+nt-1)%nt
    return xp, tp


def gather(field, index, parity, dest):
    global neighbor, lattice
    i = None
    j = None

    global EVEN, ODD, EVENANDODD, gen_pt

    if parity == EVEN:
        idx = np.where(lattice.parity == EVEN)[0]
        dest[idx] = lattice[field][neighbor[index,idx]]
        

    if parity == ODD:
        idx = np.where(lattice.parity == ODD)[0]
        dest[idx] = lattice[field][neighbor[index,idx]]
                       
    if parity == EVENANDODD:
        gen_pt = lattice[field][neighbor[index,:]].copy()


def hmc():
    global lattice, volume, nf, DELTAMAX, con, MINUS, EVENANDODD, counter

    i, j, m, n = 0, 0, 0, 0
    hold, hnew, deltah = 0.0, 0.0, 0.0
    xx = 1.0
    z = 0
    max_count = 100
    count = 0

    while count < max_count and xx>z:
        lattice.phi[:] = lattice.sigma[:]
        lattice.mom[:] = np.array(list(gasdev() for i in range(volume)))
        lattice.eta[:,:] = np.random.randn(volume,nf)

        for n in range(nf):
            matp2p("eta","chi",MINUS,EVENANDODD,n)

        hold = hamil(0,1) # initial value of the Hamiltonian
        piup(step/2) # initial half step

        # Leap frog loop for n-1 full steps
        for j in tqdm(range(0,mdstep-1)):
            lattice.phi[:] += step*lattice.mom[:]
            piup(step)

        lattice.phi[:] += step*lattice.mom[:]

        piup(step/2)

        # calculation of new Hamiltonian
        hnew = hamil(1,1)

        # Accept/reject step
        import random
        xx = random.random()
        deltah = hnew - hold
        print("The new hamil is", hnew)
        if deltah > DELTAMAX:
            print(f"HMC loop {count} REJECTED in CALL {counter} for LARGE DELTAH.\n")
            print("\n Program terminated.\n")
            for i in range(volume):
                print(f"{con[i]}\n")
            for m in range(meas_loop):
                print(f"{ac_store[m]}\n")

            system.exit(1)
        else:
            z = math.exp(-deltah)

        count += 1
    if xx <= z:
        lattice.sigma[:] = lattice.phi[:]

    return count


def layout():
    global no_odd_sites, no_even_sites, volume

    no_even_sites = volume/2
    no_odd_sites = volume/2


def site_index(x,t):
    i = None
    xr = None
    tr = None

    global nx, nt

    xr = x%nx
    tr = t%nt

    i = xr + nx*tr

    return i


def hamil(hflag, flag):
    global lattice, volume, g, nf, cgiter1, residue1
    i,n = 0,0
    h = 0

    h = np.sum(((0.5/(g*g)) * lattice.phi* lattice.phi) + (0.5 * lattice.mom * lattice.mom))

    if hflag != 0:
        for n in range(nf):
            cg_md("chi","eta",cgiter1, residue1, flag, n)
            h = np.sum(lattice.chi[:,n] * lattice.eta[:,n])
    else:
        for n in range(nf):
            h = np.sum(lattice.eta[:,n] * lattice.eta[:,n])

    return h


def zerolat():

    global lattice, volume

    print("ZEROLAT: All zero initial config. of `sigma' field.")
    lattice.sigma[:] = 0


def coldlat():

    global lattice, volume

    print(" COLDLAT: Cold initial config. of `sigma' field")

    lattice.sigma[:] = 1


def coldlat2():

    global lattice, volume

    print(" COLDLAT.4: Cold initial config. of `sigma' field")
    
    lattice.sigma[:] = 0.4


def hotlat():

    global lattice, volume

    print(" HOTLAT: Hot initial config. of `sigma' field.")

    for i in range(volume):
        lattice.sigma[i] = 2 * ran2() - 1


def filelat():

    global lattice, volume

    print(" Configuration to be read from the file sigma.in")
    print("currently not supported was feeling tired")
    # for i in range(volume):
        # lattice[i].sigma = 0.4

def funnylat():

    global lattice, volume

    for i in range(volume):
        lattice[i].sigma = i


def piup(t):
    global lattice, volume, g, nf, cgiter2, residue2, PLUS, EVENANDODD
    # print("Piup Calculation. \n")

    lattice.mom[:] -= ((1/(g*g)) * lattice.phi[:] * t)

    for n in range(nf):
        cg_md("chi","eta",cgiter2, residue2, 1, n)
        matp2d("eta","p",PLUS,EVENANDODD,n)
        lattice.mom[:] += (2 * lattice.p[:] * lattice.eta[:,n] * t)



def ran2():
    return random.random()

def gauss():
    return np.random.randn()

def gasdev():
    global iset, gset

    if iset == 0:
        while True:
            v1 = (2.0 * ran2()) - 1.0
            v2 = (2.0 * ran2()) - 1.0
            rsq = (v1 * v1) + (v2 * v2)
            if not (rsq >= 1.0 or rsq == 0.0):
                break

        fac = math.sqrt(-2.0 * math.log(rsq) / rsq)
        gset = v1 * fac
        iset = 1
        return v2 * fac
    else:
        iset = 0
        return gset

def get_lattice():
    global lattice
    return lattice
def fet_neighbour():
    global neighbor
    return neighbor

def main():
    global lattice, sw_flag, garbage, bin_length, no_garbage, ac_store, store
    #------------------------------------------------
    # File IO Defined in the code
    #------------------------------------------------
    pin = open("hmc_python\sigma1.in","r")
    ptout = open("hmc_python\sigma1.out", "a")
    ptacl = open("hmc_python\sigma1.acl", "a")
    ptlat = open("hmc_python\sigma1.lat", "a")
    ptprop = open("hmc_python\sigma1.prop", "w")
    ptpropacl = open("hmc_python\sigma1.propacl", "w")

    prompt = setup_gn()
    readin(prompt)


    # randomize() # Not defined here but present in the original code

    # filelat(pin) # Getting converted will be merged later
    coldlat() # <-----------------------------------------------Remove later !!!

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

            if no_acc > no_garbage:
                ac_store[acc] = average_sigma()
                acc += 1
                no_auto += 1

                if no_auto % seg_length == 0:
                    lbd = acc - no_auto
                    for m in range(lbd, acc):
                        seg_av_sigma += ac_store[m]/seg_length
                    autocorel(seg_av_sigma,lbd,a)
                    a += 1
                    seg_av_sigma = 0
                    no_auto = 0

            con[:] = conf[:]
            conf[:] = lattice.sigma[:]

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
                print(average_sigma())
                store[j] += average_sigma()
                t_ex_sigma += store[j]
                j += 1

            if meas%prop_length == 0:
                pass

        av_sigma = t_ex_sigma / no_meas
        sqdev = 0
        for k in range(no_meas):
            sqdev += (store[k]-av_sigma)**2
        d_av_sigma = math.sqrt(sqdev/(no_meas-1))

        acc_rate = no_acc/no_hmc

        print(f"\n\n\ no_acc_traj = {no_acc} \t no_hmc = {no_hmc} \t acc_rate = {acc_rate} \n\n")
        print(f"av_sigma={av_sigma} \t av_psi{av_psi}\n\n")
        print(f"d_av_sigma={d_av_sigma} \t d_av_psi{d_av_psi}\n\n")

        for i in range(0,volume):
            ptlat.write(f"{lattice[i].sigma}\n")

    else:
        print("Fault !!!")
