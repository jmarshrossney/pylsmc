"""
Script.

Suggests boundary positions given a number of subdomains, range, and number of bins per
subdomain (assumed to be constant).
Suggestion based on optimal bin width being ~ 2:1:2 going from low_boundary:0:high_boundary.

    python suggest_boundaries.py Ns low_boundary high_boundary Nbins_per_subdomain
"""

import numpy as np
from sys import argv

#################
## User inputs ##
#################

Ns = int(argv[1]) # Number of subdomains

bmin = abs(float(argv[2])) # Minimum mu (boundaries[0])
bmax = abs(float(argv[3])) # Maximum mu (boundaries[Ns+1])

# Assuming you've got the same number of bins per subdom...
Nbins = float(argv[4]) # Proposed number of bins per subdom


##########################################################################
## Suggest boundaries based on factor of ~2 difference in subdom widths ##
##########################################################################

offset = 0.015 * (bmax + bmin) * np.sign(bmax-bmin)

xlow = np.linspace(1,2,Ns/2)
xhigh = np.linspace(1,2,Ns/2)

Nlow = (bmin + offset) / np.cumsum(xlow)[-1]
Nhigh = (bmax - offset) / np.cumsum(xhigh)[-1]

blow = np.cumsum(xlow)*Nlow*-1 + offset
bhigh = np.cumsum(xhigh)*Nhigh + offset

blow_copy = np.copy(blow)
for i in range(Ns/2):
    blow_copy[i] = blow[Ns/2-1-i]

boundaries = np.concatenate( (blow_copy, [offset], bhigh) )

widths = np.zeros(Ns)
bin_widths = np.zeros(Ns)
for i in range(Ns):
    widths[i] = boundaries[i+1] - boundaries[i]
    bin_widths[i] = widths[i] / Nbins

ratios = widths / np.min(widths)

print "Suggested boundaries:" 
print "( " + "%1.2f, "*(Ns+1) %tuple(boundaries) + ")"

print ""

print "Subdom          Boundaries           Ratio of widths          Bin width "
for i in range(Ns):
    print "%d           %06.2f,   %06.2f            %6.2f                 %06.4f" \
            %(i,boundaries[i],boundaries[i+1],ratios[i],bin_widths[i])

