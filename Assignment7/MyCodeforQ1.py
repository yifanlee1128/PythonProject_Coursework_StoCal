import numpy                as np      # numerical computing library
import matplotlib.pyplot    as plt
import numpy.random      as ran

# initial parameters
X0=0
t=50.
dt=0.001

# ensure the number of steps is an integer
numebrOfSteps=t/dt
numebrOfSteps = int(numebrOfSteps)+1
dt = t/numebrOfSteps

recordOfXt=[]

# for loop to calculate the sample value

Xt=X0
recordOfXt.append(Xt)
for j in range(numebrOfSteps):
    Xt= Xt - Xt * dt + ran.normal(0, np.sqrt(dt), 1)[0]
    recordOfXt.append(Xt)

# calculate sample mean and standard deviation
plt.plot([x*dt for x in range(numebrOfSteps+1)],recordOfXt)

plt.xlabel("t")
plt.ylabel("Xt")
str_para='Ornstein Uhlenbeck process with delta_t='+'%.5f' % dt
plt.title(str_para)
plt.savefig('Ornstein Uhlenbeck process.pdf')
plt.show()
