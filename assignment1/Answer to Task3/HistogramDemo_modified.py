#  Python 3.7
#  For the course Stochastic Calculus, Fall semester 2019, Courant Institute, NYU
#  Author: course instructor, Jonathan Goodman
#  https://www.math.nyu.edu/faculty/goodman/teaching/StochCalc2019/StochCalc.html

#  filename: HistogramDemo.py
#  For assignment 1.

#  Create and plot a histogram of independent Gaussian random variables

import numpy             as np      # numerical computing library
import matplotlib.pyplot as plt     # library of plot routines
import numpy.random      as ran     # random number generators

print("Python file HistogramDemo.py")
print("How to make and display a histogram")


m = 1000   # number of random samples
n = 17      # dimension
L = 50      # number of bins
a = 0       # left end of the histogram region
b = 8      # right end of the histogram region



#    Set up the bins, the bin centers, and the bin count array

dx = (b-a)/L   # dx = bin size = (length being binned)/(number of bins)

N  = np.zeros( shape = [L], dtype =  int )  # bin counts, initialized to zero

#    xc = array of bin centers, uniformly spaced (using linspace).
#         The smallest one is a*.5*dx (half a bin size away from a)

xc = np.linspace( start = a+.5*dx, stop = b-.5*dx, num = L)

i=0
r=[0.*m]
while(i<n):
    tempdata=ran.randn(m)
    r+=tempdata**2
    i=i+1
r=np.sqrt(r)

for X in r:                          # look at each X in the X array
    k = int(np.floor( ( X-a )/dx ) )  # the bin X goes into
    if ( ( k >= 0 ) and ( k < L ) ):  # check that k is in the bin range
        N[k] += 1                     # increment the count for bin k
uh = N/(m*dx)                         # uh = u hat = estimate of PDF

Neb = np.ndarray( [L], dtype=float)   # N (bin count) error bar
for k in range(L):
    p = dx*uh[k]                      # estimated probability of a sample in B[k]
    Neb[k] = np.sqrt(m*p*(1.-p))      # square root of the Bernoulli variance
ueb = Neb/(m*dx)                      # standard dev of u is s.dev of N, scaled

#     Evaluate the true PDF using the standard normal density formula

rMax = 10. + n        # integrate to rMax instead of infinity.
                      # write 10. instead of 10 so it will be floating point
nPts = 1000 + 10*n*n  # number of points for integration, use a lot
dr   = rMax/(nPts-1)  # width of an integration cell.
                      # it's (nPts-1) because if nPts = 3 then there are 2 cells
sum = 0.
for r in ( np.linspace(0, rMax, nPts)):
   sum += np.power(r,(n-1))*np.exp(-r*r/2)
Z = dr*sum
pdf = lambda x: 1/Z*np.power(x,(n-1))*np.exp(-x*x/2)   # Gaussian pdf formula (1/sqrt(2*pi))e^{-x^2/2}
ue  = np.ndarray([L], dtype = float)  # ue is for "u exact" (the actual pdf)
for k in range(L):
    ue[k] = pdf(xc[k])                # evaluate the density formula at bin centers


plt.errorbar( xc, uh, yerr=ueb, linestyle = '', marker = 'o', label = 'histogram')

plt.plot( xc, ue, 'b--', label = 'R')
title = "R density estimation, " + str(m) + " samples, " + str(L) + " bins,"+str(n)+"demensions"
plt.title(title)
plt.xlabel('r')
plt.ylabel('u')

plt.legend()
plt.savefig('HistogramDemo.pdf')
plt.show()
