#  Python 3.7
#  For the course Stochastic Calculus, Fall semester 2019, Courant Institute, NYU
#  Author: course instructor, Jonathan Goodman
#  https://www.math.nyu.edu/faculty/goodman/teaching/StochCalc2019/StochCalc.html

#  filename: RandomWalkDemo.py
#  For assignment 2.

#  Simulate random walk paths up to N steps
#  Plot a histogram of hitting times
#  Compare to the ``theoretical'' density of hitting for Brownian motion

import numpy             as np  # numerical computing library
import matplotlib.pyplot as plt  # library of plot routines
import numpy.random      as ran  # random number generators

print("Assignment 2")
print("Random walk hitting times")

n_max = 3000  # maximum number of random walk steps
m = 20000  # number of paths
pu = .49  # probability that X -> X+dx in one step
pd = 1. - pu  # probability that X -> X+dx in one step
dx = .104239  # step size for stepping right
x0 = 2.  # starting point

a=(pu-pd)/dx # parameter in the condition with drift

Nk = np.zeros(n_max, dtype=int)  # Nk[k] will count times K=k

cpe = np.zeros(n_max, dtype=float)  # empirical cumulative probability
# cpe[k] = empirical prob(K<k)

cpt = np.zeros(n_max, dtype=float)  # theoretical cumulative probability
# cpt[k] = theoretical prob(K<k)
# Brownian motion approximation

#   Generate m paths and record the hitting time histogram

for path in range(m):
    X = x0
    for k in range(n_max):
        U = ran.uniform()
        if (U < pu):  # Pr( U < pu ) = pu, if U in [0,1] is uniform
            X = X + dx  # move to the right by dx
        else:  # move left with probability 1-pu = pd
            X = X - dx
        if (X < 0.):  # did the path hit zero?
            Nk[k] += 1  # record the hitting time is k
            break  # start the next path

#    Compute the CDF from the hitting times

cpe_now = 0.  # cumulative probability to the present k
cpt_now = 0.  # cpe_now = empirical, cpt_now = theoretical

dt = dx * dx  # for Brownian motion approx, dt for one step
tk = .5 * dt  # center of the first time interval

for k in range(n_max):
    cpe_now += Nk[k] / m  # divide by number of paths to get a probability
    cpe[k] = cpe_now  # record the empirical cumulative probabilities

    #  Brownian motion theoretical approximation
    # the formula has been changed

    pdf = (np.exp(-a*x0)*x0 / np.sqrt(2 * np.pi * tk * tk * tk)) * np.exp(-(x0 * x0+a*a*tk*tk) / (2 * tk))
    cpt_now += pdf * dt
    cpt[k] = cpt_now
    tk += dt

fig = plt.figure(figsize=(6, 6))  # Make one plot file with two plots

ax1 = fig.add_axes([.1, .6, 0.85, 0.35])  # locations by trial and error
ax1.plot(Nk, label='empirical')
ax1.set_title("hitting time histogram")
ax1.set_xlabel("k")
ax1.set_ylabel("N")
ax1.grid()

ax2 = fig.add_axes([.1, .1, 0.85, 0.35])
ax2.plot(cpe, label='empirical')
ax2.plot(cpt, label='theoretical')
ax2.set_title("cumulative probability")
ax2.set_xlabel("k")
ax2.set_ylabel("P")
ax2.legend()
ax2.grid()

#    put the run parameters into blank space in the CDF plot
#    put two parameters per line place them not on the grid lines

param_string = "x0 = {0:7.3f}, dx = {1:9.6f}".format(x0, dx)
ax2.text(n_max / 2., .5, param_string)

param_string = "pu = {0:7.3f}, m = {1:7.2e}".format(pu, m)
ax2.text(n_max / 2., .3, param_string)

param_string = "a= {0:9.6f}".format(a)
ax2.text(n_max / 2., .1, param_string)

plt.savefig('RandomWalkDemo_withDrift.pdf')
plt.show()
