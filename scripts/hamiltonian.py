import numpy as np
from scipy import linalg as la
import copy
import matplotlib.pyplot as plt

def vii_integrand():
    return

def tii_integrand():

def comp_simpson(left, right, num_subint, func):
    summation = 0
    dx = (float(right) - float(left)) / float(num_subint)
    i = 0
    while i < num_subint + 1:
        if i == 0 or i == num_subint:
            summation = summation + func(left + i * dx)
        if i % 2 == 0 and i != 0 and i != num_subint:
            summation = summation + 2 * func(left + i * dx)
        if i % 2 == 1 and i != num_subint:
            summation = summation + 4 * func(left + i * dx)
        i = i + 1
    area = (dx / 3) * summation
    return area

def bubble_sort(eig_energies, eig_vectors):
    new_list = copy.copy(eig_energies)
    new_mat = copy.deepcopy(eig_vectors)
    num_pairs = len(new_list) - 1
    for j in range(num_pairs):
        for i in range(num_pairs - j):
            if new_list[i] > new_list[i+1]:
                new_list[i], new_list[i+1] = new_list[i+1], new_list[i]
                new_mat[:,[i, i+1]] = new_mat[:,[i+1,i]]
    return new_list, new_mat

r = [float(x.split()[1]) for x in open('potential.out').readlines()]
potent = [float(x.split()[2]) for x in open('potential.out').readlines()]
n = len(r)
# Left boundary
r1 = r[0]
# Right boundary
r2 = r[n-1]
# reduced mass
mu = 0.5

# Create potential matrix
pe_matrix = np.zeros((n, n))
for i in range(n):
    pe_matrix[i,i] =

# Create kinetic matrix
ke_matrix = np.zeros((n, n))
ke_matrix[0][0] = -2
ke_matrix[1][0] = 1
ke_matrix[n - 1][n - 1] = -2
ke_matrix[n - 2][n - 1] = 1
for i in range(1, n - 1):
    ke_matrix[i - 1][i] = 1
    ke_matrix[i][i] = -2
    ke_matrix[i + 1][i] = 1
ke_matrix = (-1/(2*mu)) * ke_matrix

hamiltonian = ke_matrix + pe_matrix
eig_val, eig_vec = la.eig(hamiltonian)
sort_eigval, sort_eigvec = bubble_sort(eig_val, eig_vec)
f = open('2d_eigenvalues.dat', 'w+')
f.write("n       E_n\n")
for i in range(n):
    f.write("%s %s\n" % (i+1, sort_eigval[i]))
f.close()

# groundprob = np.array([sort_eigvec[i][0]*sort_eigvec[i][0] for i in range(n)])
# first_prob = np.array([sort_eigvec[i][1]*sort_eigvec[i][1] for i in range(n)])
# sec_prob = np.array([sort_eigvec[i][2]*sort_eigvec[i][2] for i in range(n)])
# third_prob = np.array([sort_eigvec[i][3]*sort_eigvec[i][3] for i in range(n)])
# plt.plot(r, groundprob, label='n = 1', color='dodgerblue')
# plt.plot(r, first_prob, label='n = 2', color='red')
# plt.plot(r, sec_prob, label='n = 3', color='green')
# plt.plot(r, third_prob, label='n = 4', color='purple')
# plt.legend(loc='upper right', fontsize=13)
# plt.xlim(0, max(r))
# plt.ylim(min(groundprob)-0.03, max(groundprob)+0.03)
# plt.xlabel("r", fontsize=15)
# plt.ylabel(r'$|\psi_n(x)|^2$', fontsize=13)
# plt.tight_layout()
# plt.savefig("prob_density.pdf")
# plt.clf()
