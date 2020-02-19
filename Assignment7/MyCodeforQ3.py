import numpy                as np      # numerical computing library
import matplotlib.pyplot    as plt
import numpy.random      as ran

# initial parameters
X0=1
t=50.
dt=0.001
a=0
b=2
# ensure the number of steps is an integer
numebrOfSteps=t/dt
numebrOfSteps = int(numebrOfSteps)+1
dt = t/numebrOfSteps

recordOfXt=[]
recordOfLt=[]
recordOfMt=[]
# for loop to calculate the sample value

Xt=X0
Lt=0
Mt=0
recordOfXt.append(Xt)
recordOfLt.append(Lt)
recordOfMt.append(Mt)
for j in range(numebrOfSteps):
    temp=Xt +ran.normal(0, np.sqrt(dt), 1)[0]
    dLk=max((a-temp),0)
    dMk=max((temp-b),0)
    Xt=temp+dLk-dMk
    Lt=Lt+dLk
    Mt=Mt+dMk
    recordOfXt.append(Xt)
    recordOfLt.append(Lt)
    recordOfMt.append(Mt)
t_span=[i*dt for i in range(numebrOfSteps+1)]
plt.plot(t_span, a*np.ones(numebrOfSteps+1), label = 'Lower boundary')
plt.plot(t_span, b*np.ones(numebrOfSteps+1), label = 'Upper boundary')
plt.plot(t_span, recordOfMt, label = 'Mt')
plt.plot(t_span, recordOfLt, label = 'Lt')
plt.plot(t_span, recordOfXt, label = 'Xt')
plt.xlabel('t')
plt.title('reflecting Brownian motion')
plt.legend()
plt.savefig('reflecting Brownian motion.pdf')
plt.show()



avg=np.mean(recordOfXt)
Y=np.array(recordOfXt-avg)
m=len(Y)
C=[]

for j in range(0,int((2/dt)),1):
    C.append(np.sum(Y[0:(m-j)]*Y[(0+j):m])/(m-j))
t_span1=[i*dt for i in range(0,int(2/dt))]
plt.semilogy(t_span1 , C, label = 'auto-covariance')
plt.legend()
title = 'auto-covariance graph with t = ' + str(t)
plt.title(title)
plt.savefig('auto-covariance function.pdf')
plt.show()
