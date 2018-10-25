import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

from sys import argv, exit

from params import *
import domain as dom

# Check for correct number of args
if len(argv) < 3:
    print "Please provide two arguments: (1) number of subdomains, (2) input diffusivity file"
    exit(1)

# Tolerance (integral of rootD for different subdoms differ by no more than tol)
tol = 0.01
Iround = str(tol)[::-1].find('.')

# Number of subdomains as argument
Ns = int(argv[1])

# Load diffusivity from file
diffusivity = np.loadtxt(argv[2])
diffusivity = diffusivity/np.min(diffusivity)

# Add ghost points to edges of diffusivity so boundaries are included
tmp = np.zeros(len(diffusivity) + 2)
tmp[1:-1] = diffusivity
tmp[0] = diffusivity[0] + 0.5*(diffusivity[0] + diffusivity[1])
tmp[-1] = diffusivity[-1] + 0.5*(diffusivity[-1] + diffusivity[-2])
diffusivity = tmp

# Need to extrapolate mu_bins too ughh
mu_bins = dom.get_global_mu_bins()
tmp = np.zeros(len(mu_bins) + 2)
tmp[1:-1] = mu_bins
tmp[0] = mu_bins[0] + 0.5*(mu_bins[0] + mu_bins[1])
tmp[-1] = mu_bins[-1] + 0.5*(mu_bins[-1] + mu_bins[-2])
mu_bins = tmp

# Interpolate rootD to create function to integrate
rootD = interp1d( mu_bins, np.sqrt(diffusivity), kind='linear')
I_rootD = np.zeros(Ns) # integral of rootD in each subdom

# Initialise subdomain boundaries - equal spacing
trial_boundaries = np.linspace(boundaries[0], boundaries[-1], Ns+1)

# Boundaries shift at each iteration
shift_frac = 0.1    # initial fraction of subdom width to shift boundaries
Nosc = 10           # Reduce shift_frac if boundaries oscillating for Nosc iterations
reduce_frac = 0.5   # If oscillating, shift_frac <-- reduce_frac * shift_frac
updown = [ [], ] * (Ns-1)   # Tells us when it's oscillating
counter = 0         # Count up to Nosc before checking if oscillating
bround = str(float(shift_frac))[::-1].find('.')

converged = False
while converged == False:
   
    counter += 1

    # Iterate over subdomains, integrating rootD between trial boundaries
    for s in range(Ns):
        
        # Take overlaps into account
        if s > 0:
            l_olap = max(abs_olap, frac_olap*(trial_boundaries[s+1]-trial_boundaries[s]) )
        else:
            l_olap = 0
        if s < Ns-1:
            h_olap = max(abs_olap, frac_olap*(trial_boundaries[s+1]-trial_boundaries[s]) )
        else:
            h_olap = 0
        mu_l = trial_boundaries[s] - l_olap
        mu_h = trial_boundaries[s+1] + h_olap
        
        integration_result = quad(rootD, mu_l, mu_h)
        I_rootD[s] = integration_result[0]


    # Shift boundaries
    for s in range(Ns-1):
        subdom_width = trial_boundaries[s+1] - trial_boundaries[s]
        if I_rootD[s] < I_rootD[s+1]:
            trial_boundaries[s+1] += subdom_width * shift_frac
            updown[s].append(1)
        elif I_rootD[s] > I_rootD[s+1]:
            trial_boundaries[s+1] -= subdom_width * shift_frac
            updown[s].append(-1)
        else:
            updown.append(0)


    # Test to see if shift_frac needs to be lowered
    if counter > Nosc:
        oscillating = True
        for s in range(Ns-1):
            if np.sum( np.array(updown[s][-Nosc:]) ) != 0:
                oscillating = False

        if oscillating == True:
            print "boundaries: "
            print np.around(np.array(trial_boundaries), bround)
            print "integral: "
            print np.around(np.array(I_rootD), Iround)
            print "shift frac %f -> %f" %(shift_frac, reduce_frac*shift_frac)
            print ""
            shift_frac = reduce_frac*shift_frac
            bround = str(float(shift_frac))[::-1].find('.')
            updown = [ [], ] * (Ns-1)
            counter = 0


    # Test for equality
    if np.max(I_rootD) - np.min(I_rootD) < tol:
        converged = True


# End while loop


print ""
print "Finished"
print "boundaries: "
print repr(np.around(np.array(trial_boundaries), bround))




