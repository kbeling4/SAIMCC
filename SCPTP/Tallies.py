import matplotlib.pyplot as plt

# ---------- tallies -------------------------
def e_tallier( tally, bins, state, initial ):
    for i in range( state[2], len(bins) ):
        if state[1] < bins[i]:
            tally[i-1] += 1 / (initial[1] - state[1])
            state[2] = i
            break
    return tally

def p_tallier( tally, bins, state ):
    for i in range( state[2], len(bins) ):
        if state[0] < bins[i]:
            tally[i-1] += 1
            state[2] = i
            break
    return tally, state

def normalizer( tally, bins, initial ):
    for i in range( 0, len(tally) ):
        b = (bins[i + 1] + bins[i]) / 2
        tally[i] = tally[i] * ( initial[1] - b )
    return tally

def plotter( tally, bins, typ ):
    x = []
    y = []
    for i in range( 0, len(bins) - 1 ):
        x.append( (bins[i+1] + bins[i])/2 )
        y.append( tally[i] )
    if typ == 'bar':
        w = bins[2] - bins[1]
        plt.bar( x, y, width=w)
        plt.xlim( bins[0], bins[-1] )
        plt.title('1700 MeV protons on 0.5cm W')
        plt.xlabel('Emerging Energy (MeV)')
        plt.ylabel('Fraction in Energy bin')
    if typ == 'scatter':
        plt.plot( x, y, 'o')
        plt.xlim( bins[0], bins[-1] )
    plt.show()

    
        
