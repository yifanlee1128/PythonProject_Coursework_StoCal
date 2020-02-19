import numpy                as np      # numerical computing library
import matplotlib.pyplot    as plt
import numpy.random      as ran

# initial parameters
X0=0
t=20.
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

avg=np.mean(recordOfXt)
Y=np.array(recordOfXt-avg)
m=len(Y)
C=[]

for j in range(0,int((2/dt)),1):
    C.append(np.sum(Y[0:(m-j)]*Y[(0+j):m])/(m-j))

t_span=[i*dt for i in range(0,int(2/dt))]
C_1 = 1/2 * np.exp(-np.array(t_span))
plt.plot(t_span,C,label = 'auto-covariance')
plt.plot(t_span,C_1,label = 'theoretical auto-covariance')
plt.legend()
title = 'auto-covariance graph with t = ' + str(t)
plt.title(title)
plt.savefig('auto-covariance function t=20.pdf')
plt.show()


