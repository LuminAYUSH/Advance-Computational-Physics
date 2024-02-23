import hmc_python
from hmc_python.util import *
from hmc_python.lparam import *

def main():
    global lattice
    #------------------------------------------------
    # File IO Defined in the code
    #------------------------------------------------
    pin = open("C:/Users/cdipt/Documents/GitHub/Advance-Computational-Physics/hmc_python/sigma1.in","r")
    ptout = open("C:/Users/cdipt/Documents/GitHub/Advance-Computational-Physics/hmc_python/sigma1.out", "a")
    ptacl = open("C:/Users/cdipt/Documents/GitHub/Advance-Computational-Physics/hmc_python/sigma1.acl", "a");
    ptlat = open("C:/Users/cdipt/Documents/GitHub/Advance-Computational-Physics/hmc_python/sigma1.lat", "a");
    ptprop = open("C:/Users/cdipt/Documents/GitHub/Advance-Computational-Physics/hmc_python/sigma1.prop", "w");
    ptpropacl = open("C:/Users/cdipt/Documents/GitHub/Advance-Computational-Physics/hmc_python/sigma1.propacl", "w");

    prompt = setup_gn()
    readin(prompt)

    print(lattice)
    print(get_lattice())