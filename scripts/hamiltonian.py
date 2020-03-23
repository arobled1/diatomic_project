import numpy as np
from scipy import linalg as la
import copy
import matplotlib.pyplot as plt

def vij_integrand(left, right, potential, indexi, indexj, num_points, red_mass, quant_num):
    constant1 = 2/(right - left)
    deltax = (float(right) - float(left)) / float(num_points - 1)
    r = np.arange(left, right+deltax, deltax)
    centrifugal = (quant_num*(quant_num+1)) / (2*red_mass*r**2)
    sinelist_1 = np.sin ( (indexi*np.pi*(r - left) ) / (right - left) )
    sinelist_2 = np.sin ( (indexj*np.pi*(r - left) ) / (right - left) )
    return constant1 * (potential + centrifugal) * sinelist_1 * sinelist_2

def tii_integrand(left, right, index, red_mass, num_points):
    constant2 = ( (index * np.pi)**2 ) / (2 * red_mass * (right - left)**2)
    return constant2

def comp_simpson(left, right, num_subint, list):
    summation = 0.0
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

def sine_basis(left, right, coeffs, x):
    basis_sum = 0.0
    for i in range(len(coeffs)):
        basis_sum = basis_sum + coeffs[i] * np.sqrt(2/(right - left)) * np.sin ( ((i+1) * np.pi*(x - left) ) / (right - left) )
    return basis_sum

r = [float(x.split()[1]) for x in open('potential.out').readlines()]
r = np.asarray(r)
potent = [float(x.split()[2]) for x in open('potential.out').readlines()]
potent = np.asarray(potent)
subintervals = len(potent) - 1
n = len(r)
# Left boundary
r1 = r[0]
# Right boundary
r2 = r[n-1]
# reduced mass
mu = 0.957069
# vlaue of l
l = 5

# Create potential matrix
pe_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        # v_integrand = vii_integrand(r1, r2, potent, i+1, j+1, n)
        v_integrand = vij_integrand(r1, r2, potent, i+1, j+1, n, mu, l)
        pe_matrix[i,j] = comp_simpson(r1, r2, subintervals, v_integrand)

ke_matrix = np.zeros((n, n))
for i in range(n):
    ke_matrix[i,i] = tii_integrand(r1, r2, i+1, mu, n)

hamiltonian = ke_matrix + pe_matrix
eig_val, eig_vec = la.eig(hamiltonian)

sort_eigval, sort_eigvec = bubble_sort(eig_val, eig_vec)
if np.sign(sort_eigvec[0,0]) == -1:
    sort_eigvec = sort_eigvec * -1
f = open('2d_eigenvalues_l5.dat', 'w+')
f.write("n       E_n\n")
for i in range(n):
    f.write("%s %s\n" % (i+1, sort_eigval[i]))
f.close()

ground = []
for j in range(len(r)):
    ground.append(sine_basis(r1, r2, sort_eigvec[:,0], r[j]))
second = []
for j in range(len(r)):
    second.append(sine_basis(r1, r2, sort_eigvec[:,1], r[j]))
third = []
for j in range(len(r)):
    third.append(sine_basis(r1, r2, sort_eigvec[:,2], r[j]))
fourth = []
for j in range(len(r)):
    fourth.append(sine_basis(r1, r2, sort_eigvec[:,3], r[j]))
# twenty = []
# for j in range(len(r)):
    # twenty.append(sine_basis(r1, r2, sort_eigvec[:,19], r[j]))

# plt.plot(r, ground, label='n = 1', color='dodgerblue')
# plt.plot(r, second, label='n = 2', color='red')
# plt.plot(r, third, label='n = 3', color='green')
# plt.plot(r, fourth, label='n = 4', color='purple')
# # plt.plot(r, twenty, label='n = 20', color='purple')
# plt.legend(loc='upper right', fontsize=13)
# plt.xlim(0.5, 3.0)
# plt.ylim(min(ground)-1, max(ground)+1)
# plt.xlabel("r", fontsize=15)
# plt.ylabel(r'$\psi_n(r)$', fontsize=13)
# plt.tight_layout()
# plt.savefig("wavefunction_l1.pdf")
# plt.clf()

ground_prob = np.array([ground[i]*ground[i] for i in range(n)])
second_prob = np.array([second[i]*second[i] for i in range(n)])
third_prob = np.array([third[i]*third[i] for i in range(n)])
fourth_prob = np.array([fourth[i]*fourth[i] for i in range(n)])
# twenty_prob = np.array([twenty[i]*twenty[i] for i in range(n)])
plt.plot(r, ground_prob, label='n = 1', color='dodgerblue')
plt.plot(r, second_prob, label='n = 2', color='red')
plt.plot(r, third_prob, label='n = 3', color='green')
plt.plot(r, fourth_prob, label='n = 4', color='purple')
# plt.plot(r, twenty_prob, label='n = 20', color='purple')
plt.legend(loc='upper right', fontsize=13)
plt.xlim(0.5, 3.0)
plt.ylim(min(ground_prob)-1, max(ground_prob)+1)
plt.xlabel("r", fontsize=15)
plt.ylabel(r'$|\psi_n(r)|^2$', fontsize=13)
plt.tight_layout()
plt.savefig("prob_density_l5.pdf")
plt.clf()
