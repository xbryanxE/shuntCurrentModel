"""
Check that a solution satisfies boundary constraints for the corresonding boundaries
"""
def linear_u_test(uh_, R, V, bounds, idx):
    Rch = R[0]
    Rm = R[1]
    r = 0
    for i in range(2):
        
        if idx[i] > 0:
            if i == 0:
                bounds_ = (uh_.x.array[1] + V[0]) / (2 * Rch + Rm)    
            else:
                bounds_ = (uh_.x.array[-2] + V[-1]) / (2 * Rch + Rm)

            r += (bounds_ - bounds[i])**2
            bounds[i] = bounds_
    return [bounds, r]
