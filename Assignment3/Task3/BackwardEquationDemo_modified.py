#  Python 3.7
#  For the course Stochastic Calculus, Fall semester 2019, Courant Institute, NYU
#  Author: course instructor, Jonathan Goodman
#  https://www.math.nyu.edu/faculty/goodman/teaching/StochCalc2019/StochCalc.html

#  filename: BackwardEquationDemo.py
#  For assignment 3.

#  Use a finite difference method to compute the (approximate) solution
#  of the backward equation for a 1D Brownian motion with absorbing boundaries
#  at x=a and x=b.

#  Illustrate what happens if you try to run the backward equation the
#  wrong way

import numpy                as np  # numerical computing library
import matplotlib.pyplot    as plt  # library of plot routines
from matplotlib.animation import FFMpegWriter  # for making a movie

print("Assignment 3")
print("Brownian motion backward equation finite difference method")

plt.rcParams['animation.ffmpeg_path'] = 'C:/Users/liyf4/PycharmProjects/StochasticCalculus/ffmpeg/bin/ffmpeg.exe'

xl = -1.  # solve on the interval [xl,xr]
xr = 1.  # "xl" for "left" endpoint, etc.
T = 1.  # five final conditions at time T
a = -10.  # drift velocity


#     The function that computes the payout function

def v(x):
    """compute and return the payout function"""
    return x * x


n = 100  # number of points in the interior of the interval [a,b]
lam_d = .4  # desired CFL ratio. The ratio used may be slightly smaller

dx = (xr - xl) / (n + 1)  # distance between grid points
xLocs = np.linspace(xl + dx, xr - dx, n)

dt_d = 2 * lam_d * dx * dx  # time step using desired value of lambda
nt_d = T / dt_d  # number of time steps with the desired lambda,
# probably not an integer
nt = int(nt_d) + 1  # round up the number of time steps, round down dt
dt = T / nt  # the time step to use
lam = .5 * dt / (dx * dx)  # corresponding value of lambda for the computation

#             weights for the finite difference time step

pp = lam+a*dt/2/dx
p0 = 1. - 2 * lam
pm = lam-a*dt/2/dx

#             solution arrays and initialization with final condition

fk = np.ndarray(n)  # the solution at time t_k
fkp1 = np.ndarray(n)  # the solution at time t_{k+1}

for j in range(n):
    x = xl + (j + 1) * dx
    fk[j] = v(x)

#     Set up to record the movie

MovieSeconds = 10  # number of seconds the movie should take
fps = int(nt / MovieSeconds)  # frames/sec = (frames)/(seconds)

metadata = dict(title='Backward equation dynamics', artist='Matplotlib',
                comment='Demo')
writer = FFMpegWriter(fps=fps, metadata=metadata)

fig = plt.figure()  # dunno what this does, but it works for me
l, = plt.plot([], [])

plt.xlim(xl, xr)  # you have to set the same limits for each frame
plt.ylim(0., 1.)
plt.grid()

with writer.saving(fig, "BackwardEquationMovie_modified.mp4", 100):
    #    The main finite difference time step loop
    tk = T
    for k in range(nt):  # nt is the number of time steps

        #        The first and last points are special
        #        The absorbing boundary condition gives f = 0 at a and b

        fkp1[0] = p0 * fk[0] + pp * fk[1]  # the first point, nobody on the left
        fkp1[n - 1] = pm * fk[n - 2] + p0 * fk[n - 1]  # the last point, nobody on the right

        #        This is the main finite difference formula

        for j in range(1, (n - 1)):  # the interior points go from 1 to n-2
            fkp1[j] = pm * fk[j - 1] + p0 * fk[j] + pp * fk[j + 1]

        #         The newly computed values will be the old values in the next iteration

        fk[:] = fkp1[:]  # fk = fkm1 does the wrong thing.
        tk -= dt

        l.set_data(xLocs, fk)  # plot the solution at time tk
        title = "Value function at time {0:7.3f}, a = {1:6.2f}".format(tk, a)
        plt.title(title)
        plt.xlabel("x")
        plt.ylabel("f")
        writer.grab_frame()  # make this plot a movie frame

