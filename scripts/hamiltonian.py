import numpy as np
from scipy import linalg as la
import copy
import matplotlib.pyplot as plt

def vii_integrand(left, right, potential, index, num_points):
    constant1 = 2/(right - left)
    deltax = (float(right) - float(left)) / float(num_points - 1)
    r = np.arange(left, right+deltax, deltax)
    sinelist_v = (index*np.pi*(r - left) ) / (right - left)
    sinelist_v = sinelist_v**2
    return constant1 * potential * sinelist_v

def tii_integrand(left, right, index, red_mass, num_points):
    constant2 = ( (index * np.pi)**2 ) / (red_mass * (right - left))
    deltax = (float(right) - float(left)) / float(num_points - 1)
    r = np.arange(left, right+deltax, deltax)
    sinelist_t = (index*np.pi*(r - left) ) / (right - left)
    sinelist_t = sinelist_t**2
    return constant2 * sinelist_t

def comp_simpson(left, right, num_subint, list):
    summation = 0
    dx = (float(right) - float(left)) / float(num_subint)
    i = 0
    while i < num_subint + 1:
        if i == 0 or i == num_subint:
            summation = summation + list[i]
        if i % 2 == 0 and i != 0 and i != num_subint:
            summation = summation + 2 * list[i]
        if i % 2 == 1 and i != num_subint:
            summation = summation + 4 * list[i]
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
potent = np.asarray(potent)
subintervals = len(potent) - 1
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
    v_integrand = vii_integrand(r1, r2, potent, i, n)
    pe_matrix[i,i] = comp_simpson(r1, r2, subintervals, v_integrand)

t_matrix = np.zeros((n, n))
for i in range(n):
    t_integrand = tii_integrand(r1, r2, i, mu, n)
    t_matrix[i,i] = comp_simpson(r1, r2, subintervals, t_integrand)

# hamiltonian = ke_matrix + pe_matrix
# eig_val, eig_vec = la.eig(hamiltonian)
# sort_eigval, sort_eigvec = bubble_sort(eig_val, eig_vec)
# f = open('2d_eigenvalues.dat', 'w+')
# f.write("n       E_n\n")
# for i in range(n):
#     f.write("%s %s\n" % (i+1, sort_eigval[i]))
# f.close()

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
