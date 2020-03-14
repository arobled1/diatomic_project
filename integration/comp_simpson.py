import numpy as np
import matplotlib.pyplot as plt

#===============================================================================
# Uses composite simpson's rule to approximate the integral of a sin(x)^2 from
# [a,b].
# a = 1 and b = 10 in this code.
# This code runs composite simpson's rule iteratively where each iteration adds
# one more subinterval until the stopping criteria is met.
# The stopping criteria is set so the code stops when the absolute value of the
# difference between the approximate area of the current and previous iteration
# is less than some epsilon.
# The true value of the integral (using Mathematica) is 4.499088044.
#
# By: Alan Robledo
# Updated date: March 14, 2020
#===============================================================================
# Output:
#       area = 4.502341569241636
#       Number of iterations: 15
#===============================================================================

def func(xin):
    return np.sin(xin)**2

def comp_simpson(a, b):
    num_subint = 1
    area_new = 0
    area_prev = 0
    j = 0
    conv = []
    iter = []
    while j == 0:
        area_prev = area_new
        area_new = 0
        summation = 0
        dx = (float(b) - float(a)) / float(num_subint)
        i = 0
        while i < num_subint + 1:
            if i == 0 or i == num_subint:
                summation = summation + func(a + i * dx)
            if i % 2 == 0 and i != 0 and i != num_subint:
                summation = summation + 2 * func(a + i * dx)
            if i % 2 == 1 and i != num_subint:
                summation = summation + 4 * func(a + i * dx)
            i = i + 1
        area_new = (dx / 3) * summation
        conv.append(abs(area_new - area_prev))
        iter.append(num_subint)
        if abs(area_new - area_prev) <= epsilon:
            j = 1
        num_subint = num_subint + 1
    print("\nApproximate area under the curve is: %s" % area_new)
    print("Precision = |area_new - area_previous| = %s" % abs(area_new - area_prev))
    print("Numer of iterations required: %s" % iter[len(iter) - 1])
    return iter, conv

a = 1
b = 10
epsilon = 10**-2
iter, conv = comp_simpson(a, b)

plt.xlim(min(iter) - 1, max(iter) + 1)
plt.ylim(min(conv) - 0.1, max(conv) + 0.1)
plt.plot(iter, conv, '-o')
plt.xlabel("# of iterations")
plt.ylabel("|area_current - area_previous|")
plt.savefig("integral_convergence.pdf")
plt.clf()
print("\nCheck your directory for graphical verification of convergence!")
