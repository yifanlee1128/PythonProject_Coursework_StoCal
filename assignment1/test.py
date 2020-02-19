
import numpy             as np      # numerical computing library
import matplotlib.pyplot as plt     # library of plot routines
import numpy.random      as ran     # random number generators
import math
m=10000

Xa = ran.randn(m)    # an array of n independent standard normals.
print(math.sqrt(np.sum(np.square(Xa))))
R=[math.sqrt(np.sum(np.square(Xa))) for i in range(m)]
print(R)