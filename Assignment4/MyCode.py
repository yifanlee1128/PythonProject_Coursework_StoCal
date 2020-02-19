import numpy                as np      # numerical computing library
import matplotlib.pyplot    as plt
import numpy.random      as ran

# initial parameters
N = 1000
S0=1
sigma=1.0
r=0.0
t=1.
dt=0.005
L=30

# ensure the number of steps is an integer
numebrOfSteps=t/dt
numebrOfSteps = int(numebrOfSteps)+1
dt = t/numebrOfSteps

recordOfST=[]

# for loop to calculate the sample value
for i in range(N):
    St=S0
    for j in range(numebrOfSteps):
        St=St+r*St*dt+sigma*St*ran.normal(0,np.sqrt(dt),1)[0]
    recordOfST.append(St)

# calculate sample mean and standard deviation
mean_St=np.mean(recordOfST)
std_St=np.std(recordOfST)
print("The sample mean is: ",mean_St)
print("The sample standard deviation is: ",std_St)

# plot histogram
recordOfST.sort()
a= recordOfST[0]-0.5 #-10
b= recordOfST[-1]+0.5 #20
dx = (b-a)/L
Nm  = np.zeros( shape = [L], dtype =  int )
for X in recordOfST:                          # look at each X in the X array
    k = int(np.floor( ( X-a )/dx ) )  # the bin X goes into
    if ( ( k >= 0 ) and ( k < L ) ):  # check that k is in the bin range
        Nm[k] += 1
uh = Nm/(N*dx)
xc = np.linspace( start = a+.5*dx, stop = b-.5*dx, num = L)
title="Histogram od St: s0="+str(S0)+",sigma="+str(sigma)+",r="+str(r)+",t="+str(t)
plt.bar( xc, uh, label = 'histogram')
plt.title(title)
plt.xlabel("St")
plt.ylabel("density")
plt.show()