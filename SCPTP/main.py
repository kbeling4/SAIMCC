# Simple 1D monte carlo charged particle transport code with energy in MeV and position in cm
import numpy as np
import random
import time
from threading import Thread

import Loop as loop
import Tallies as tal
# --------- loop function ----------------------------

def Main():
    t0 = time.time()
    # num_threads = 2
    nps = 1e6
    random.seed( 123456 )

    # --------- problem properties ----------------
    # Physical properties
    start  =  0   # Units of cm
    thickness = 0.1

    initialEnergy  = 1700 # Units of MeV
    energyCutOff = 1e-1

    particle = ['proton', 1,  938.27231]
    material = ['Tungsten', 74, 183.84, 19.3, 6e-04]

    # Tally properties
    num_bins = 100
    bins = np.linspace( 1690, 1700, num=num_bins)
    tally = [0]*(num_bins - 1)

    # --------- run problem --------------------

    initial = [ start, initialEnergy, thickness, energyCutOff ]

    loop.SCPT(particle, material, bins, tally, nps, initial)

    # threads = [None]*int(num_threads)
    # for i in range(0, int(num_threads)):
    #     nps_t = nps/num_threads
    #     t = Thread( target=loop.SCPT, args=(particle, material, bins, tally, nps_t, initial) )
    #     t.start()
    #     threads[i] = t
        
    # for t in threads:
    #     t.join()
    
    # -------- results -------------------------
    
    print "nps = %s" % int(nps)
    t1 = time.time()
    print '---------------------'
    total = t1 - t0
    print "total time = %s" % total

    tal.normalizer( tally, bins, initial, nps )
    tal.plotter( tally, bins, 'scatter' )


# ---------- Run problem ---------------------

if __name__ == '__main__':
    Main()
