import numpy as np

import Constants as const
import Newton as newt

def Qmax( T, particle ):
    Gamma = ( T + particle[2] ) / particle[2]
    Beta = 1 - ( 1 / Gamma )
    qmax = 1.022 * Beta * Gamma
    return [ Beta, Gamma, qmax ]

def TotXS( T, particle, material ):
    Amp = 0.1536 * particle[1]**2 * material[1] * material[3] / material[2]
    BGQ = Qmax( T, particle )
    return Amp * ( 1 / BGQ[0] ) * ( ( 1 / material[4] - 1 / BGQ[2] ) - ( BGQ[0] / BGQ[2] )
                                        * np.log( BGQ[2] / material[4] ) )

def Pusher( state, rand, particle, material ):
    XS = TotXS( state[1], particle, material )
    delx = ( 1 / XS ) * np.log(1/rand)
    state[0] += delx
    return state

def Decrementer( state, rand, particle, material ):
    q = Qmax( state[1], particle )
    delq = ( ( ( 1 - rand ) / material[4]) + ( rand / q[2] ) )**(-1)
    state[1] -= delq
    return state

def dfdQ( state, particle, material, Q ):
    qm = Qmax( state[1], particle )
    XS = ( 1/material[4] - 1/qm[2] ) - (qm[0]/qm[2])*np.log(qm[2]/material[4])
    return ( 1/(Q*Q) - (qm[0]/qm[2])*(1/Q) )/XS

def f_Q( state, particle, material, rand, Q ):
    qm = Qmax( state[1], particle )
    XS = ( 1/material[4] - 1/qm[2] ) - (qm[0]/qm[2])*np.log(qm[2]/material[4])
    return (( 1/material[4] - 1/Q ) - (qm[0]/qm[2])*np.log(Q/material[4]))/XS - rand

def Dec_E( state, rand, particle, material ):
    delQ = newt.Newton( f_Q, dfdQ, 1e-3, rand, 1e-3, state, particle, material )
    state[1] -= delQ
    return state
 
