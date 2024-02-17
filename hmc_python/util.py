### Util for the HMC code everything is non-vectorised in this Util
import numpy as np
import scipy
import sklearn as sk
import matplotlib.pyplot as plt
import sys
import os
import math


class psferm:
    def __init__(self):
        self.f = None # Has to be initialized as a array of size [4]


class site:
    def __init__(self):
        self.x = None
        self.t = None
        self.parity = None
        self.sigma = None
        self.phi = None
        self.mom = None
        self.chi = psferm()
        self.eta = psferm()
        self.p = None
        self.r = None
        self.mp = None
        self.mmp = None

    def __getitem__(self,key):
        return self.__dict__[key]
    def __setitem__(self,key,value):
        self.__dict__[key] = value


def setup_gn():
    prompt = initial_set()
    layout();
    make_lattice();
    make_nn_gather();
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

    status = getprompt() #Using the direct value inplace of the pointer in the original code base

    if status !=0:
        print("error in input: initial prompt")
        return -1

    nx = get_i(prompt, "nx")
    nt = get_i(prompt, "nt")

    if nx%2 !=0 and nt%2 !=0:
        print("nx, nt must be even!! \n")
        sys.exit()

    print(f"Lattice dimensions = {nx} {nt}\n")

#     Switch flag
    sw_flag = get_i(prompt,"switch_flag")
    print(f"Switch_Flag = {sw_flag}\n")

#     Number of measurements
    no_garbage = get_i(prompt,"no_of_garbage_loops")
    print(f"No of garbage loops = {no_garbage}")

    # the length of each bin
    bin_length = get_i(prompt,"bin_length")
    print(f"bin_length = {bin_length}")

    # Number of HMC iterations, only accepted ones count
    hmc_it = get_i(prompt,"no_of_hmc_iterations")
    print(f"# of hmc iterations = {hmc_it}")

    # length after which fermionic observables are measured
    meas_length = get_i(prompt,"meas_length")
    print(f"meas_length = {meas_length}")

    # length after which the propagator autocorrelations are calculated
    prop_length = get_i(prompt,"prop_length")
    print(f"prop_length = {prop_length}")

    # segment length after which autocorrelations are measured
    seg_length = get_i(prompt,"seg_length")
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
    global gen_pt
    global T_int
    global T_int_prop

    i = None
    j = None
    k = None

    x = None
    t = None

    no_bin = no_garbage/bin_length
    no_meas = meas_loop/meas_length
    no_a_seg = meas_loop/seg_length
    no_prop_seg = meas_loop/prop_length

    #All allocations managed by python backend
    lattice = [site()]*volume
    store = [None]*no_meas
    conf = [None]*volume
    con = [None]*volume
    garbage = [None]*no_garbage
    ac_store = [None]*meas_loop
    ac_prop = [None]*prop_length
    bin_av = [None]*no_bin
    psi = [None]*no_meas
    psi_acl = [None]*meas_loop
    G_prop = [None]*nt
    G_store = [[None]*no_meas]*nt # This should be nx see the defination using LSIZE
    G_temp = [[None]*meas_loop]*nt  # This should be nx see the defination using LSIZE
    prop = [None]*nt
    tprop = [None]*nt

    gen_pt = [None]*volume

    T_int = [[[None]*no_a_seg]*MAXT_cut]*NOT_cut
    T_int_prop = [[[[None]*no_prop_seg]*nt]*MAXT_cut]*NOT_cut

    for t in range(nt):
        for x in range(nx):
            i = site_index(x,t)     # Function not defined yet !!!!
            lattice[i].x = x
            lattice[i].t = t
            if t%2 ==0:
                lattice[i].sign = 1
            else:
                lattice[i].sign = -1
                if (x+t)%2 == 0:
                    lattice[i].parity = EVEN;
                else:
                    lattice[i].parity = ODD

    for i in range(meas_loop):
        ac_store[i] = 0
    for i in range(no_bin):
        bin_av[i] = 0
    for i in range(nt):
        G_prop = 0
    for i in range(NOT_cut):
        for j in range(MAXT_cut):
            for k in range(no_a_seg):
                T_int[i][j][k] = 0


def get_f(prompt,variable_name_string):

    if prompt == 1:
        x = double(input(f"enter {variable_name_string}"))
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

    g = get_f(prompt,"g")
    nf = get_i(prompt,"nf")

    mdstep = get_i(prompt,"no_of_md_steps")
    step = get_f(prompt,"step_size")
    print(f"no of md steps = {mdstep}, step size = {step}")

    cgiter1 = get_i(prompt,"max_cg_iterations_for_hamil")
    cgiter2 = get_i(prompt,"max_cg_iterations_for_piup")
    printf(f"maximum no. of conj. grad. iterations = {cgiter1},{cgiter2}")

    residue1=get_f(prompt,"residue_for_cg_hamil");
    residue2=get_f(prompt,"residue_for_cg_piup");
    printf(f"residues for conjugate grad = {residue1},{residue2}")


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


def average_sigma():
    global volume
    i = None
    av_sigma = None
    t_sigma = None
    t_sigma = 0

    for i in range(volume):
        t_sigma += lattice[i].sigma
    av_sigma = t_sigma/volume
    return av_sigma


def F_PT(site,fo):
    global lattice
    "Here site is the index of the global variable lattice"
    return list(lattice[site].__dict__.keys())[fo]


def F_OFFSET(a):
    return list(site().__dict__.keys()).index(a)


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

    for i in range(volume):
        dsize_src += (lattice[i][F_PT(i,src)].f[flavor])**2

    size_src = math.sqrt(dsize_src)

    if cgflag == 0:
        for i in range(volume):
            lattice[i][F_PT(i,dest)].f[flavor] = 0
            lattice[i].r = lattice[i][F_PT(i,src)].f[flavor]
            lattice[i].p = lattice[i].r

        dsize_r = 1
        size_r = dsize_r

    if cgflag != 0:
        matp2d(dest,"p",PLUS,EVENANDODD,flavor)
        matd2d("p","mp",MINUS,EVENANDODD)

        dsize_r = 0

        for i in range(volume):
            lattice[i].r = lattice[i][F_PT(i,src)].f[flavor] - lattice[i].mp
            dsize_r += (lattice[i].r)**2
            lattice[i].p = lattice[i].r

        size_r = math.sqrt(dsize_r/size_src)

    cp = dsize_r

    N_iter = 0

    while N_iter < cgiter and size_r > residue:
        c = cp

        matd2d("p","mp",PLUS,EVENANDODD)

        d=0

        for i in range(volume):
            d += (lattice[i].mp)**2

        a = c/d

        matd2d("mp","mmp",MINUS,EVENANDODD)

        cp = 0
        for i in range(volume):
            lattice[i][F_PT(i,dest)].f[flavor] += a*lattice[i].p
            lattice[i].r -= a*lattice[i].mmp
            cp += (lattice.r)**2

        b = cp/c
        dsize_r = 0
        for i in range(volume):
            lattice[i].p = lattice[i].r + b*lattice[i].p
            dsize_r += (lattice[i].r)**2
        size_r = math.sqrt(dsize_r)/size_src

        N_iter = N_iter + 1


    if size_r > residue:
        print("CG_MD Not Converged")
        system.exit(1)


def cg_prop(src , dest , cgiter , residue, cgflag):
    global volume, lattice, PLUS, MINUS, EVENANDODD
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
    for i in range(volume):
        dsize_src += lattice[i][F_PT(i,src)]**2
    size_src = math.sqrt(dsize_src)

    # Initial guess
    if cgflag == 0:
        for i in range(volume):
            lattice[i][F_PT(i,dest)] = 0.0
            lattice[i].r = lattice[i][F_PT(i,src)]
        dsize_r = 1.0
        size_r = dsize_r

    if cgflag != 0:
        matd2d(dest,"mp",PLUS,EVENANDODD)
        dsize_r = 0.0
        for i in range(volume):
            lattice[i].r = lattice[i][F_PT(i,src)] - lattice[i].mp
            dsize_r += (lattice[i].r)**2
        size_r = math.sqrt(dsize_r)/size_src

    matd2d("r","p",MINUS,EVENANDODD)

    cp = 0.0
    for i in range(volume):
        cp += (lattice[i].p)**2

    # Start of CG iteration loop
    while N_iter < cgiter and size_r > residue:
        c = cp

        matd2d("p","mp",PLUS,EVENANDODD)

        d = 0.0
        for i in range(volume):
            d += (lattice[i].mp)**2

        a = c/d

        for i in range(volume):
            lattice[i][F_PT(i,dest)] += a*lattice[i].p
            lattice[i].r -= a*lattice[i].mp

        matd2d("r","mp",MINUS,EVENANDODD)

        cp = 0.0
        for i in range(volume):
            cp += (lattice[i].mp)**2

        b = cp/c

        dsize_r = 0.0
        for i in range(volume):
            lattice[i].p = lattice[i].mp + b*lattice[i].p
            dsize_r += (lattice[i].r)**2

        size_r = math.sqrt(dsize_r)/size_src
        N_iter += 1

    if size_r > residue:
        print("CG_PROP Not Converged")


def matd2d(src,dest,isign,parity):
    i = None
    n = None

    global XUP,XDN,TUP,TDN, EVENANDODD, gen_pt, volume

    gather(src, XUP, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)] = lattice[i].sign * 0.5 * gen_pt[i]

    gather(src, XDN, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)] -= lattice[i].sign * 0.5 * gen_pt[i]

    gather(src, TUP, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)] +=  0.5 * gen_pt[i]

    gather(src, TDN, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)] -=  0.5 * gen_pt[i]

    for i in range(volume):
        lattice[i][F_PT(i,dest)] = isign*(lattice[i][F_PT(i,dest)]) + (lattice[i].phi)*(lattice[i][F_PT(i,src)])


def matd2p(src,dest,isign,parity,flavor):
    i = None
    n = None

    global volume, gen_pt, XUP, XDN, TDN, TUP, EVENANDODD

    gather(src, XUP, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)].f[flavor] = lattice[i].sign * 0.5 * gen_pt[i]

    gather(src, XDN, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)].f[flavor] -= lattice[i].sign * 0.5 * gen_pt[i]

    gather(src, TUP, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)].f[flavor] += 0.5 * gen_pt[i]

    gather(src, TDN, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)].f[flavor] -= 0.5 * gen_pt[i]

    for i in range(volume):
        lattice[i][F_PT(i,dest)].f[flavor] = isign*lattice[i][F_PT(i,dest)].f[flavor] + lattice[i].phi*lattice[i][F_PT(i,src)]


def matp2d(src,dest,isign,parity,flavor):
    i = None
    n = None

    global volume, gen_pt, XUP, XDN, TDN, TUP, EVENANDODD

    gather(src, XUP, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)] = lattice[i].sign * 0.5 * gen_pt[i].f[flavor]

    gather(src, XDN, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)] -= lattice[i].sign * 0.5 * gen_pt[i].f[flavor]

    gather(src, TUP, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)] += 0.5 * gen_pt[i].f[flavor]

    gather(src, TDN, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)] -= 0.5 * gen_pt[i].f[flavor]

    for i in range(volume):
        lattice[i][F_PT(i,dest)] = isign*lattice[i][F_PT(i,dest)] +lattice[i].phi*lattice[i][F_PT(i,src)].f[flavor]


def matp2p(src,dest,isign,parity,flavor):
    i = None
    n = None

    global XUP,XDN,TUP,TDN, EVENANDODD, gen_pt, volume

    gather(src, XUP, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)].f[flavor] = lattice[i].sign * 0.5 * gen_pt[i].f[flavor]

    gather(src, XDN, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)].f[flavor] -= lattice[i].sign * 0.5 * gen_pt[i].f[flavor]

    gather(src, TUP, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)].f[flavor] +=  0.5 * gen_pt[i].f[flavor]

    gather(src, TDN, EVENANDODD, gen_pt)
    for i in range(volume):
        lattice[i][F_PT(i,dest)].f[flavor] -=  0.5 * gen_pt[i].f[flavor]

    for i in range(volume):
        lattice[i][F_PT(i,dest)].f[flavor] = isign*(lattice[i][F_PT(i,dest)].f[flavor]) + (lattice[i].phi)*(lattice[i][F_PT(i,src)].f[flavor])


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

    tprop = [0.0] * nt

    for i in range(volume):
        lattice[i].r = 0.0

    for t in range(nt):
        source = (t * nx) + mid
        lattice[source].r = 1.0

        cg_prop("r", "mmp", cgiter1, residue1, 0)

        for x in range(nx):
            xl = x * nx
            xu = ((x + 1) * nx) - 1
            n0 = x - t
            if n0 < 0:
                n0 += nt
            prop[n0] = 0.0

            for m in range(xl, xu):
                prop[n0] += lattice[m].mmp

            tprop[n0] += prop[n0]

        lpsi += lattice[source].mmp

    for t in range(nt):
        tprop[t] = tprop[t] / nt

    lpsi = lpsi / nt

    return lpsi

neighbor = [[None]*4]*volume # neighbour stores the index of the neighbouring site in place of the pointer


def make_nn_gather():
    global neighbor
    i = None
    j = None
    dir = None

    xpt = None
    tpt = None
    mul = None

    for dir in range(0,4):
        for i in range(volume):
            xpt,tpt = neighbor_coords(lattice[i].x,lattice[i].t, dir)
            j = site_index(xpt,tpt)
            neighbor[dir][i] = j


def neighbor_coords(x,t,dir):
    global neighbor
    xp = None
    tp = None

    global nx, nt

    if dir == 0:
        xp = (x+1)%nx
    if dir == 3:
        xp = (x+nx-1)%nx
    if dir == 1:
        tp = (t+1)%nt
    if dir == 2:
        tp = (t+nt+1)%nt
    return xp, tp


def gather(field, index, parity, dest):
    global neighbor
    i = None
    j = None

    global EVEN, ODD, EVENANDODD

    if parity == EVEN:
        for j in range(volume):
            if lattice[j].parity == EVEN:
                dest[j] = lattice[neighbor[index][j]][field]

    if parity == ODD:
        for j in range(volume):
            if lattice[j].parity == odd:
                dest[j] = lattice[neighbor[index][j]][field]

    if parity == EVENANDODD:
        for j in range(volume):
                dest[j] = lattice[neighbor[index][j]][field]


def hmc():
    global lattice, volume, nf, DELTAMAX, con

    i, j, m, n = 0, 0, 0, 0
    hold, hnew, deltah = 0.0, 0.0, 0.0
    xx = 1.0
    z = 0
    max_count = 100
    count = 0

    while count < max_count and xx>z:
        for i in range(volume):
            lattice[i].phi = lattice[i].sigma
            lattice[i].mom = gasdev()
            for n in range(nf):
                lattice[i].eta.f[n] = gauss()

        hold = hamil(0,1) # initial value of the Hamiltonian
        piup(step/2) # initial half step

        # Leap frog loop for n-1 full steps
        for j in range(0,mdstep-1):
            for i in range(volume):
                lattice[i].phi += step*lattice[i].mom
                piup(step)
        for i in range(volume):
            lattice[i].phi += step*lattice[i].mom
        piup(step/2)

        # calculation of new Hamiltonian
        hnew = hamil(1,1)

        # Accept/reject step
        import random
        xx = random.random()
        deltah = hnew - hold

        if deltah > DELTAMAX:
            print(f"HMC loop {count} REJECTED in CALL {counter} for LARGE DELTAH.\n")
            print("\n Program terminated.\n")
            for i in range(volume):
                print(f"{con[i]}\n")
            for m in range(meas_loop):
                print(f"{ac_store[m]}\n")

            terminate(1)
        else:
            z = math.exp(-deltah)

        count += 1

    if xx <= z:
        for i in range(volume):
            lattice[i].sigma = lattice[i].phi

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

    i = xr + nx + tr

    return i


def hamil(hflag, flag):
    global lattice, volume, g, nf, cgiter1, residue1
    i,n = 0,0
    h = 0

    for i in range(volume):
        h += (((0.5/(g*g)) * lattice[i].phi* lattice[i].phi) + (0.5 * lattice[i].mom * lattice[i].mom))

    if hflag != 0:
        for n in range(nf):
            cg_md("chi","eta",cgiter1, residue1, flag, n)
            for i in range(volume):
                h += lattice[i].chi.f[n] * lattice[i].eta.f[n]
    else:
        for n in range(nf):
            for i in range(volume):
                h += lattice[i].eta.f[n] * lattice[i].eta.f[n]

    return h


def zerolat():
    print("")
