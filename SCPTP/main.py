# Simple 1D monte carlo charged particle transport code with energy in MeV and position in cm
import numpy as np
import random
import matplotlib.pyplot as plt
import time
from threading import Thread

import Loop as loop
# --------- loop function ----------------------------

def Main():
    t0 = time.time()
    num_threads = 2
    
    # --------- problem properties ----------------
    
    nps = 1e4
    x0  =  0   # Units of cm
    e0  = 1700 # Units of MeV

    x_end = 0.1
    e_cut = 1e-1

    particle = ['proton', 1,  938.27231]
    material = ['Tungsten', 74, 183.84, 19.3, 6e-04]
    initial = [ x0, e0, x_end, e_cut ]

    num_bins = 100
    bins = np.linspace( 1690, 1700, num=num_bins)
    tally = [0]*num_bins
    random.seed( 123456 )

    num_collisions = [None]*int(nps)

    # --------- run problem --------------------

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
    
    print "nps = %s" % nps
    t1 = time.time()
    print '---------------------'
    total = t1 - t0
    print "total time = %s" % total

    # print sum(tally)
    plt.plot( bins, tally, 'o'  )
    plt.show()

# ---------- Run problem ---------------------

if __name__ == '__main__':
    Main()
