import numpy as np
import matplotlib.pyplot as plt

from params import *
import domain as dom

if algorithm == 'wang_landau': import wang_landau as alg
elif algorithm == 'multicanonical': import multicanonical as alg
elif algorithm == 'transition': import transition_matrix as alg

# Number of data points per subdomain on scatter plot
points_per_sub = 10000

if TRAP == True:
    Ns_eff = Ns
else:
    Ns_eff = 1

# ----------- #
#  Make plot  #
# ----------- #
plt.rc('text', usetex=True)
font = {'family' : 'serif',
        'size' : 14}
plt.rc('font', **font)

fig, ax = plt.subplots()
ax.set_title("Active lattice energy")
ax.set_xlabel(r"$\mu \times \beta/N$")
ax.set_ylabel(r"$E \times \beta/N$")

# Colours distinguish between subdomains
colours = ('tab:gray', 'k')

# Initialise arrays for averaging energy
E_mean = np.zeros(np.sum(bins))
counts = np.zeros(np.sum(bins))

# ----------------------------------- #
#  Iterate over different subdomains  #
# ----------------------------------- #
for s in range(Ns_eff):

    # Load the series data
    input_file = alg.file_names('pcomb', s)['s']
    input_data = np.loadtxt(input_file)

    # Unpack
    subdom = np.array(input_data[:,0], dtype=np.int)
    mu = input_data[:,1]
    weights = input_data[:,2]
    E_active = input_data[:,3]

    # Decide how frequently we want to add a point to the scatter plot
    Nrows = len(subdom)
    addpoint = int(Nrows / points_per_sub)

    # Initialise lists for scatter data
    E_scatter = []
    mu_scatter = []

    # ----------------------- #
    #  Iterate over all rows  #
    # ----------------------- #
    for row in range(Nrows):

        # Skip rows which separate independent simulations
        if subdom[row] == -1: continue

        # Ignoring overlaps, subdomain for this mu might change
        s_no_olap = dom.get_subdomain(mu[row])
        
        # Add to array using global index
        iglob = dom.get_global_index(mu[row], s_no_olap)
        E_mean[iglob] += E_active[row]
        counts[iglob] += 1

        # Add raw data to list for scatter plot
        if row % addpoint == 0:
            E_scatter.append(E_active[row])
            mu_scatter.append(mu[row])
            

    # -------------------------- #
    #  Add scatter data to plot  #
    # -------------------------- #
    mu_scatter = np.array(mu_scatter) * B / float(Natoms)
    E_scatter = np.array(E_scatter) * B / float(Natoms)
    ax.scatter(mu_scatter, E_scatter, s=1., c=colours[s%2])


# ------------------ #
#  Add mean to plot  #
# ------------------ #
mu_bins = dom.get_global_mu_bins() * B / float(Natoms)
E_mean = (E_mean / counts) * B / float(Natoms)
ax.plot(mu_bins, E_mean, 'r-')

print "Mu = E(%s) - E(%s)" %(alpha_type, beta_type)

plt.show()


