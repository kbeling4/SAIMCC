import numpy as np

def Newton(f, dfdx, x, rand, eps, state, particle, material):
    
    f_value = f( state, particle, material, rand, x)
    iteration_counter = 0
    while abs(f_value) > eps and iteration_counter < 100:
        try:
            x = x - float(f_value)/dfdx(state, particle, material, x )
        except ZeroDivisionError:
            print "Error! - derivative zero for x = ", x
            sys.exit(1)

        f_value = f( state, particle, material, rand, x)
        iteration_counter += 1

    if abs(f_value) > eps:
        iteration_counter = -1
    return x
