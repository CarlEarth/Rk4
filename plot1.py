import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

Y = np.loadtxt('time_position')
x=Y.T

ax = plt.axis([0,50,-1,1])
l=plt.plot(x[0],x[1])
plt.xlabel('t')
plt.ylabel('position')
plt.show()
