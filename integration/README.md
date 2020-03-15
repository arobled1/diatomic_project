## Numerical integration using Simpson's method

The code in 'comp_simpson.py' shows how to numerically integrate any given function (in this case, sin(x)^2) from some interval [a,b]. The code specifically uses composite simpson's rule which requires that n be an even integer in order to obtain meaningful results.

The code also shows how many subintervals are required to meet a convergence criteria. The criteria in this case being that the absolute error of the result for a given number of subintervals be less than or equal to 10^-5. The code also generates a figure that shows the absolute error as a function of the number of even subintervals.

![int_conv](https://github.com/arobled1/diatomic_project/blob/master/integration/integral_convergence.pdf)
