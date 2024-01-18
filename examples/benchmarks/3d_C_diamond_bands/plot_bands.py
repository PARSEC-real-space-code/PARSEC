#!/usr/bin/env python3

#
# Usage: python plot_bands.py
#

import matplotlib.pyplot as plt
import re
import numpy as np


with open('bands.dat', 'r') as ifs:
    data = ifs.read()

# ===========================================
# Calculate the number of k-points
#
# We count how many lines of the k-points
# Each k-point is listed as
#     2     1     1     0.00000000     0.00000000     0.00000000
# ===========================================
kps = re.findall(3*'\s+(\d+)'+3*'\s+(\S+)'+'\n', data)

print(kps[0])

num_of_kps = len(kps)

# ===========================================
# Get the eigenvalues at each k-point
#
# In bands.dat, each k-point and its eigenvalues are printed out
# in the order as specified in parsec.in
#
# Eigenvalues at each k-point are listed line by line as
# (i-th eigenvalue)    (eigenvalue in Ry)   (eigenvalue in eV)
#
# For example,
#      5         -0.627939         -8.543618
# ===========================================
eigs = re.findall('\s+(\d+)\s+(\S+)\s+(\S+)\n', data)

print(eigs[0])

num_of_eigs_per_kp = len(eigs)//num_of_kps

# Change to a numpy array (otherwise they are strings)
eigs = np.float_(eigs)

# ===========================================
# Plot the band structure in eV
# ===========================================
energy_shift = 1.911603
rydberg_in_ev = 13.605698066
x = [k for k in range(1,num_of_kps+1) for kk in range(num_of_eigs_per_kp)]
plt.scatter(x, (eigs[:,1]-energy_shift)*rydberg_in_ev, s=1)
plt.show()
