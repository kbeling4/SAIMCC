# ---------- tallies -------------------------
def e_tallier( tally, bins, state ):
    for i in range( state[2], len(bins) ):
        if state[1] < bins[i]:
            tally[i-1] += 1
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

