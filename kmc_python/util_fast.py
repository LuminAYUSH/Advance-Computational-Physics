### Util for the KMC code with numpy vectorisation :>
import numpy as np
import sys
import math
import random
from tqdm.auto import tqdm
import json

# Init File for the KMC Python package
Xup = 0
Tup = 1
Tdn = 2
Xdn = 3
Ndirs = 4

COLDLAT = 0
HOTLAT = 1
PLUS = 1
MINUS = -1

FRESH = 0
FORGET = 1
FRESHHOT = 2
MAXFILENAME = 256

lattice = None

Nx = None
Nt = None
N = None
Volume = None
Nf = None
G = None
Step = None
Eps = None
KRwarmup = None
KRmax = None
KRtraj = None
CgiterH = None
CgiterP = None
ResidueH = None
ResidueP = None

sites_on_node = None

readflag = None
saveflag = None
read_binary_sigma = None # Although defined as char array we take string
save_binary_sigma = None # Although defined as char array we take string

Nmeas = None
Measlen = None
av_sigma = None #Defined as pointer check later !!!!

iseed = None
junk_id = None

N_POINTERS = 4

gen_pt = None
neighbor = None # neighbour stores the index of the neighbouring site in place of the pointer

#-----------------------
IM1 = 2147483563
IM2 = 2147483399
AM = 1.0 / IM1
IMM1 = IM1 - 1
IA1 = 40014
IA2 = 40692
IQ1 = 53668
IQ2 = 52774
IR1 = 12211
IR2 = 3791
NTAB = 32
NDIV = 1 + IMM1 / NTAB
EPS = 1.2e-7
RNMX = 1.0 - EPS
#-----------------------
# These are the static variable defined inside the codebase
iset = 0
gset = None
iseed = -772
iseed2 = 123456789
iy = 0
iv = [0] * NTAB


class LATTICE:
    def __init__(self,volume_):
        self.x = np.array(list(0 for i in range(volume_)))
        self.t = np.array(list(0 for i in range(volume_)))
        self.index = np.array(list(0 for i in range(volume_)))
        self.sign = np.array(list(0 for i in range(volume_)))

        self.mom = np.array(list(0.0 for i in range(volume_)))
        self.pi = np.array(list(0.0 for i in range(volume_)))
        self.sigma = np.array(list(0.0 for i in range(volume_)))
        self.phi = np.array(list(0.0 for i in range(volume_)))
        self.rho = np.array(list(0.0 for i in range(volume_)))

        self.eta = np.array(list(list(0.0 for i in range(4)) for j in range(volume_)))
        self.zeta = np.array(list(list(0.0 for i in range(4)) for j in range(volume_)))

        self.omega = np.array(list(list(0.0 for i in range(4)) for j in range(volume_)))
        self.xi = np.array(list(list(0.0 for i in range(4)) for j in range(volume_)))

        self.p = np.array(list(0.0 for i in range(volume_)))
        self.r = np.array(list(0.0 for i in range(volume_)))
        self.mp = np.array(list(0.0 for i in range(volume_)))
        self.mmp = np.array(list(0.0 for i in range(volume_)))

    def __getitem__(self,key):
        return self.__dict__[key]
    def __setitem__(self,key,value):
        self.__dict__[key] = value

class params:
    def __init__(self):
        self.stopflag = None
        self.iseed = None
        self.nx = None
        self.nt = None
        self.nflv = None
        self.g = None
        self.step = None
        self.eps = None
        self.krmax = None
        self.krtraj = None
        self.cgiter1 = None
        self.cgiter2 = None

        self.residue1 = None
        self.residue2 = None

        self.read_sigma = None
        self.save_sigma = None

        self.read_bin_double = None
        self.save_binary_double = None

par_buf = params()


def initial_set():

    prompt = getprompt()

    # print("GN model with HMC algorithm\n")
    # print("Alpha machine, Version 1\n")
    # print("type 0 for reading from JSON, 1 for prompts\n")

    prompt = getprompt() #Using the direct value inplace of the pointer in the original code base

    return prompt

def readin(prompt,params_dict = None):
    global Nf, G, Step, Eps, KRtraj, KRmax, Measlen, CgiterH, ResidueH, CgiterP, ResidueP, KRwarmup
    global Nx, Nt, N, Volume
    
    if prompt == 1:
        Nf = get_i(prompt, "num_flavors")
        G = get_lf(prompt, "coupling")
        Step = get_lf(prompt, "step_size")
        Eps = get_lf(prompt, "noise_parameter")
        KRtraj = get_i(prompt, "kr_trajectory_length")
        KRmax = get_i(prompt, "kr_refreshing_steps")
        KRwarmup = get_i(prompt, "kr_warmup_steps")
        Measlen = get_i(prompt, "measurement_interval")
        CgiterH = get_i(prompt, "max_cg_iter_for_hamil")
        ResidueH = get_lf(prompt, "residue_cg_hamil")
        CgiterP = get_i(prompt, "max_cg_iter_for_piup")
        ResidueP = get_lf(prompt, "residue_cg_piup")
        Nx = get_i(prompt, "nx")
        Nt = get_i(prompt, "nt")
    
    elif prompt == 0 or prompt == 2:
        if params_dict:
            Nf = params_dict["num_flavors"]
            G = params_dict["coupling"]
            Step = params_dict["step_size"]
            Eps = params_dict["noise_parameter"]
            KRtraj = params_dict["kr_trajectory_length"]
            KRmax = params_dict["kr_refreshing_steps"]
            KRwarmup = params_dict["kr_warmup_steps"]
            Measlen = params_dict["measurement_interval"]
            CgiterH = params_dict["max_cg_iter_for_hamil"]
            ResidueH = params_dict["residue_cg_hamil"]
            CgiterP = params_dict["max_cg_iter_for_piup"]
            ResidueP = params_dict["residue_cg_piup"]
            Nx = params_dict["nx"]
            Nt = params_dict["nt"]
        else:
            print("ERROR ACCESSING THE JSON file")
    else:
        print("Invalid Prompt")
    
    if Nx == Nt:
        N = Nx
    else:
        print("Nx and Nt must be same, exiting")
        sys.exit()

    Volume = Nx*Nt
    # print(f"Lattice dimensions = {Nx} {Nt}\n")

def make_lattice():

    global lattice, Volume, Nx, neighbor

    lattice = LATTICE(Volume)
    neighbor = np.array(list(list(0.0 for i in range(Volume)) for j in range(4)))

    for t_ in range(Nt):
        for x_ in range(Nx):
            i = site_index(x_,t_)     # Function not defined yet !!!!
            lattice.x[i] = x_
            lattice.t[i] = t_
            lattice.index[i] = x_+Nx*t_
            if t_%2 ==0:
                lattice.sign[i] = 1
            else:
                lattice.sign[i] = -1

def setup_layout():
    global Volume, sites_on_node
    sites_on_node = Volume

def get_lf(prompt,variable_name_string):

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
#     prompt = int(input("Enter the prompt type (only 1 is supported for now):"))
    prompt = 1
    return prompt

def average_sigma():
    global Volume, lattice

    av_sigma = None
    t_sigma = None
    t_sigma = 0

    t_sigma = np.sum(lattice.sigma[:])
    av_sigma = t_sigma/Volume
    return av_sigma

    
def average_pbp():

    print("check this function with sir")
    
    
def congrad(src,dest,cgiter,residue,cgflag,flv):
    niter = None
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

    global Volume
    global lattice

    dsize_src = 0

    dsize_src = np.sum((lattice[src][:, flv])**2)
    size_src = math.sqrt(dsize_src)

    if cgflag == 0:

        lattice[dest][:,flv] = 0
        lattice.r[:] = lattice[src][:,flv]
        lattice.p[:] = lattice.r[:]

        dsize_r = 1
        size_r = dsize_r

    if cgflag != 0:
        matp2d(dest, "mp", PLUS, flv)
        matd2d("mp", "mmp", MINUS)

        dsize_r = 0

        lattice.r[:] = lattice[src][:,flv] - lattice.mmp[:]
        dsize_r = np.sum(lattice.r[:]**2)
        lattice.p[:] = lattice.r[:]

        size_r = math.sqrt(dsize_r)/size_src

    cp = dsize_r

    niter = 0

    while niter < cgiter and size_r > residue:
        c = cp

        matd2d("p","mp",PLUS)
        matd2d("mp","mmp",MINUS)

        d=0
        d = np.sum((lattice.p[:]*lattice.mmp[:]))
        a = c/d

        cp = 0
        lattice[dest][:,flv] = lattice[dest][:,flv] + a*lattice.p[:]
        lattice.r[:] = lattice.r[:] - a*lattice.mmp[:]
        cp = np.sum(lattice.r[:]**2)

        b = cp/c
        dsize_r = 0

        lattice.p[:] = lattice.r[:] + b*lattice.p[:]
        dsize_r = np.sum(lattice.r[i]**2)
        size_r = math.sqrt(dsize_r)/size_src
        # print(f"CG_MD Step no {niter}",size_r)

        niter = niter + 1

    if size_r > residue:
        print("CG_MD Not Converged")
        print(size_r)
        system.exit(1)

def matd2d(src,dest,isign):
    i = None
    n = None

    global Xup,Xdn,Tup,Tdn, volume, lattice

    gather(src, Xup, gen_pt)
    lattice[dest][:] = 0.5*lattice.sign[:]*gen_pt[:]

    gather(src, Xdn, gen_pt)
    lattice[dest][:] = lattice[dest][:] - 0.5*lattice.sign[:] * gen_pt[:]

    gather(src, Tup, gen_pt)
    lattice[dest][:] = lattice[dest][:] + 0.5 * gen_pt[:]

    gather(src, Tdn, gen_pt)
    lattice[dest][:] = lattice[dest][:] - 0.5 * gen_pt[:]

    lattice[dest][:] = isign*(lattice[dest][:]) + (lattice.phi[:])*(lattice[src][:])


def matd2p(src, dest, isign, flavor):
    i = None
    n = None

    global Xup, Xdn, gen_pt, Tdn, Tup, lattice

    gather(src, Xup, gen_pt)
    lattice[dest][i,flavor] = lattice.sign[:] * 0.5 * gen_pt[:]

    gather(src, Xdn, gen_pt)
    lattice[dest][:,flavor] = lattice[dest][:,flavor] - 0.5*lattice.sign[:] * gen_pt[:]

    gather(src, Tup, gen_pt)
    lattice[dest][:,flavor] = lattice[dest][:,flavor] + 0.5 * gen_pt[:]

    gather(src, Tdn, gen_pt)
    lattice[dest][:,flavor] = lattice[dest][:,flavor] - 0.5 * gen_pt[:]

    lattice[dest][:,flavor] = isign*lattice[dest][:,flavor] + lattice.phi[:]*lattice[src][:]


def matp2d(src,dest,isign,flavor):

    global Xup,gen_pt, Xdn, Tdn, Tup, lattice

    gather(src, Xup, gen_pt)
    lattice[dest][:] = 0.5*lattice.sign[:] * gen_pt[:,flavor]

    gather(src, Xdn, gen_pt)
    lattice[dest][:] = lattice[dest][:] - 0.5 * lattice.sign[:] * gen_pt[:,flavor]

    gather(src, Tup, gen_pt)
    lattice[dest][:] =  lattice[dest][:] + 0.5 * gen_pt[:,flavor]

    gather(src, Tdn, gen_pt)
    lattice[dest][:] = lattice[dest][:] - 0.5 * gen_pt[:,flavor]

    lattice[dest][:] = isign*lattice[dest][:] +lattice.phi[:]*lattice[src][:,flavor]



def matp2p(src,dest,isign,flav):

    global Xup, Xdn, Tup, Tdn, gen_pt, lattice

    gather(src, Xup, gen_pt)
    lattice[dest][:,flav] = lattice.sign[:] * 0.5 * gen_pt[:,flav]

    gather(src, Xdn, gen_pt)
    lattice[dest][:,flav] -= lattice.sign[:] * 0.5 * gen_pt[:,flav]

    gather(src, Tup, gen_pt)
    lattice[dest][:,flav] +=  0.5 * gen_pt[:,flav]

    gather(src, Tdn, gen_pt)
    lattice[dest][:,flav] -=  0.5 * gen_pt[:,flav]


    lattice[dest][:,flav] = isign*(lattice[dest][:,flav]) + (lattice.phi[:])*(lattice[src][:,flav])





def nn_coords(x,t,dir):
    global neighbor
    xp = x
    tp = t

    global Nx, Nt

    if dir == 0:
        xp = (x+1)%Nx
    if dir == 3:
        xp = (x+Nx-1)%Nx
    if dir == 1:
        tp = (t+1)%Nt
    if dir == 2:
        tp = (t+Nt-1)%Nt
    return xp, tp

def make_nn_gathers():

    #This will run only once hence not vectorised

    global neighbor, lattice, Volume
    i = None
    j = None
    dir = None

    xpt = None
    tpt = None
    mul = None

    for dir in range(0,4):
        for i in range(Volume):
            xpt, tpt = nn_coords(lattice.x[i],lattice.t[i], dir)
            j = site_index(xpt,tpt)
            neighbor[dir,i] = j
    neighbor = neighbor.astype(np.uint16)


def gather(field, index, dest):

    global neighbor, lattice
    global gen_pt
    gen_pt = lattice[field][neighbor[index,:]].copy()


def kramer():
    global lattice, Nf, MINUS, Eps, Step

    max_count = 100
    count = 0

    lattice.phi[:] = lattice.sigma[:]
    lattice.rho[:] = np.random.normal(loc = 0, scale = 1, size = lattice.rho.shape)

    term1 = -1*Eps*Step
    term2 = math.sqrt(1 - math.exp(2*term1))
    
    lattice.pi[:] = math.exp(term1)*lattice.mom[:] + term2*lattice.rho[:]

    # try:
    hold = hamil(1,1)
    # print("Old Hamiltonian",hold)
    piup(Step/2)
    # print("Initial Half Step taken")
    lattice.phi[:] = lattice.phi[:] + Step*lattice.pi[:]
    # print("Leapfrog done")
    piup(Step/2)
    # print("Final Half step done")   
    hnew = hamil(1,1)
    # print("New Hamiltonian",hnew)

    xx = random.random()
    deltah = hnew-hold
    if deltah > 50: #Possible explosion, all values discarded
        return -1
    # print("The value of deltah is",deltah)
    z = min(1,math.exp(-1*deltah))
    
    if xx<z:
        # print("Accepted!")
        # print("--------------------------")
        lattice.sigma[:] = lattice.phi[:]
        lattice.mom[:] = lattice.pi[:]
        return 1
    else:
        # print("Rejected!")
        # print("--------------------------")
        lattice.mom[:] = -1*lattice.pi[:]
        return 0
        


def layout():
    global no_odd_sites, no_even_sites, volume

    no_even_sites = volume/2
    no_odd_sites = volume/2


def site_index(x,t):
    global Nx, Nt
    i = None
    xr = None
    tr = None

    xr = x%Nx
    tr = t%Nt

    i = xr + Nx*tr

    return i


def hamil(hflag, flag):
    global lattice, Volume, G, Nf, CgiterH, ResidueH
    i,n = 0,0
    h = 0

    hml = np.sum(((0.5/(G*G)) * lattice.phi* lattice.phi)) + np.sum((0.5 * lattice.pi * lattice.pi))

    if hflag != 0:
        for n in range(Nf):
            congrad("zeta","eta",CgiterH, ResidueH, flag, n)
            hml = hml + np.sum(lattice.zeta[:,n] * lattice.eta[:,n])
    else:
        for n in range(Nf):
            hml = hml + np.sum(lattice.zeta[:,n] * lattice.eta[:,n])

    return hml


def zerolat():

    global lattice, volume

    # print("ZEROLAT: All zero initial config. of `sigma' field.")
    lattice.sigma[:] = 0


def coldlat():

    global lattice, volume

    # print(" COLDLAT: Cold initial config. of `sigma' field")

    lattice.sigma[:] = 1


def coldlat2():

    global lattice, volume

    # print(" COLDLAT.4: Cold initial config. of `sigma' field")

    lattice.sigma[:] = 0.4


def hotlat():

    global lattice, volume

    print(" HOTLAT: Hot initial config. of `sigma' field.")

    for i in range(volume):
        lattice.sigma[i] = 2 * ran2() - 1


def filelat():

    global lattice, volume
      # for i in range(volume):
        # lattice[i].sigma = 0.4

def funnylat():

    global lattice, volume

    for i in range(volume):
        lattice[i].sigma = i

def piup(t):
    global lattice, Volume, G, Nf, CgiterP,ResidueP, PLUS
    # print("Piup Calculation. \n")

    lattice.pi[:] = lattice.pi[:] - ((1/(G*G)) * lattice.phi[:] * t)

    for n in range(Nf):
        congrad("zeta","xi", CgiterP, ResidueP, 1, n)
        matp2p("xi","omega",PLUS,n)
        lattice.pi[:] = lattice.pi[:] + (2 * lattice.omega[:, n] * lattice.xi[:, n] * t)


        
        
def main(prompt = None, file = None, dict_=None, suppress_output = False, write_file = None):
    """
    Entry module for the code
    prompt (int): The way the parameters will be entered, 0 will read from JSON 1 will take the input from USER, 2 will take a dict
    file (str): Path to the JSON file if the parameters will be read
    write_file (str): Name of the numpy array in which the measurements will stored for analysis
    """
    global KRtraj,KRwarmup, lattice, Measlen
    console = sys.stdout
    
    configurations = []
    
    if suppress_output is False:
        print(" ===================== KRAMER'S ALGORITHM =====================")
        print(" =============Shifting from console to out.dat file============")

    AVERAGE_SIGMA = []
    AVERAGE_XI = []
    params_dict = None
    
    if prompt is None:
        prompt = initial_set()
    if prompt == 0: 
        if file is None:
            print("Path to JSON file expected")
            print("Exiting code...")
            return 1
        else:
            file = open(file, 'r')
            params_dict = json.load(file)
    if prompt == 2:
        if dict_ is None:
            print("ERROR empty dictionary passed ...")
            print("Exiting...")
            return 1
        else:
            params_dict = dict_
            
    readin(prompt,params_dict)
    setup_layout()
    make_lattice()
    make_nn_gathers()
    
    if suppress_output is False:
        f = open("./OUTPUT.dat","w")
        sys.stdout = f

    coldlat()
    l=0
    numkr=0
    numaccp=0
    
    for j in range(0,KRtraj+KRwarmup):

        lattice.mom = np.random.normal(loc=0,scale=1,size=lattice.mom.shape)
        lattice.zeta = np.random.normal(loc=0,scale=(1/math.sqrt(2)),size=lattice.zeta.shape)

        
        k = 0
        while (k<KRmax):
            krflag = kramer()
            if krflag == 0:
                numkr = numkr + 1
            if krflag == 1:
                numkr = numkr + 1
                numaccp = numaccp + 1
                k = k + 1
            if krflag == -1:
                lattice.mom = np.random.normal(loc=0,scale=1,size=lattice.mom.shape)
                lattice.zeta = np.random.normal(loc=0,scale=(1/math.sqrt(2)),size=lattice.zeta.shape)
                
        if suppress_output is False:        
            print(f"Kritter = {j}, av_sigma = {average_sigma()}")
            
        if j > KRwarmup:
            AVERAGE_SIGMA.append(average_sigma())
            if j%Measlen == 0:
                configurations.append(lattice.sigma[:][:,None].copy())
    
    acceptance_rate = numaccp/numkr
    configurations = np.concatenate(configurations,axis = -1)
    if write_file:
        np.save(write_file+".npy",configurations)
    
    if suppress_output is False:
        print(f"Acceptance rate = {acceptance_rate}")
        print(" =======================Success============================")
        print(" =============Shifting from out.dat to console ============")
        f.close()
        sys.stdout = console
        print(f"Acceptance rate = {acceptance_rate}")
    return AVERAGE_SIGMA, acceptance_rate
