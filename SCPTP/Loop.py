import random

import Constants as const
import CollisionPhysics as XS
import Tallies as tal

def SCPT( particle, material, bins, tally, nps, initial ):
    for i in range( 0, int(nps) ):
        if i % int(nps/10) == 0:
            print "nps = %s" % i
        j = 0
        state = [ initial[0], initial[1], 0 ]
    
        while j != 1:
            if state[0] >= initial[2] or state[1] <= initial[3]:
                # num_collisions[i] = k
                tal.e_tallier( tally, bins, state, initial )
                break
            else:
                # Position Mover
                r1  = random.random()
                XS.Pusher( state, r1, particle, material )

                # Energy Decrementer
                r2  = random.random()
                XS.Dec_E( state, r2, particle, material )
