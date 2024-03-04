import math

# Init File for the HMC Python package
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

XUP = 0
TUP = 1
TDN = 2
XDN = 3
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
neighbor = None # neighbour stores the index of the neighbouring site in place of the pointer

# These are the static variable defined inside the codebase
iset = 0
gset = None
