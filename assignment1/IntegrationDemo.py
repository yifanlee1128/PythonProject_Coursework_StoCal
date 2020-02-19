#  Python 3.7
#  For the course Stochastic Calculus, Fall semester 2019, Courant Institute, NYU
#  Author: course instructor, Jonathan Goodman
#  https://www.math.nyu.edu/faculty/goodman/teaching/StochCalc2019/StochCalc.html

#  filename: IntegrationDemo.py
#  For assignment 1.

#  Compute the integral in assignment 1 numerically

import numpy             as np      # numerical computing library

print("Python file IntegrationDemo.py")
print("Compute an integral")

#   Z = integral from 0 to infinite r^{n-1}e^{-r^2/2} dr

n    = 4              # dimension of the space
rMax = 10. + n        # integrate to rMax instead of infinity.
                      # write 10. instead of 10 so it will be floating point
nPts = 1000 + 10*n*n  # number of points for integration, use a lot
dr   = rMax/(nPts-1)  # width of an integration cell.
                      # it's (nPts-1) because if nPts = 3 then there are 2 cells

sum = 0.
for r in ( np.linspace(0, rMax, nPts)):
   sum += np.power(r,(n-1))*np.exp(-r*r/2)
Z = dr*sum

output = "integral for n = {0:3d} is Z = {1:12.6e}".format(n, Z)
print(output)
