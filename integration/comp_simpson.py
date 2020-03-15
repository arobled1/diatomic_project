#===============================================================================
# Uses composite simpson's rule to numerically integrate sin(x)^2 from [a,b].
# a = 1 and b = 10 is used for this code and num_subint must be even.
# The code runs composite simpson's rule iteratively where each iteration adds
# two more subintervals until the stopping criteria is met.
# The stopping criteria is set so the code stops when the absolute error is less
# than some epsilon.
# The true value of the integral is (9/2) - [sin(20) - sin(2)]/4 , which is
# approximately equal to 4.49908804402451.
#
# By: Alan Robledo
# Updated date: March 14, 2020
#===============================================================================
# Input: a = left interval
#        b = right interval
#        true = true value of integral
#        epsilon = absolute error criteria
#===============================================================================
# Output:
#       area = 4.4990785192359555
#       Number of subintervals required: 16
#===============================================================================

import numpy as np
import matplotlib.pyplot as plt

def func(xin):
    return np.sin(xin)**2

def comp_simpson(left, right, num_subint):
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

def iterate(left, right, eps):
    true = 4.49908804402451
    num_subint = 2
    area_new = 0
    j = 0
    conv = []
    iter = []
    while j == 0:
        area_new = comp_simpson(left, right, num_subint)
        conv.append( abs(area_new - true)  )
        iter.append(num_subint)
        if abs(area_new - true) <= eps:
            j = 1
        num_subint = num_subint + 2
    print("\nApproximate area under the curve is: %s" % area_new)
    print("Absolute error = %s" % abs(area_new - true))
    print("Numer of subintervals required: %s" % iter[len(iter) - 1])
    return iter, conv

a = 1
b = 10
epsilon = 10**-5
iter, conv = iterate(a, b, epsilon)

plt.xlim(min(iter) - 1, max(iter) + 1)
plt.ylim(-0.02, 0.02)
plt.plot(iter, conv, '-o')
plt.xlabel("# of subintervals")
plt.ylabel("Absolute error")
plt.savefig("integral_convergence.pdf")
plt.clf()
print("\nCheck your directory for graphical verification of convergence!")
