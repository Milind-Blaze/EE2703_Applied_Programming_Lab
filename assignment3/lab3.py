from __future__ import division

import numpy as np
from matplotlib import pyplot as plt



#functions

def exponential(x):
    return np.exp(x)

def coscos(x):
    return np.cos(np.cos(x))

valuerange= np.arange(-2*np.pi, 4*np.pi, 0.1)

plt.figure(1)
plt.grid(True)
plt.plot(valuerange, exponential(valuerange),"ro")
plt.xlabel("$x$")
plt.ylabel("$e^x$")
plt.show()
