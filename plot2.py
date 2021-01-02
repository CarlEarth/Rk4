import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches

Y = np.loadtxt('time_position')
x=Y.T
fig, ax = plt.subplots()
ax = plt.axis([0,50,-1,1])
l, =plt.plot(x[0][0:1],x[1][0:1])
#redDot, = plt.plot(x[0][0], x[1][0], 'ro')
print(len(x[0]))
def animate(i):
    #l=plt.plot(x[0][0:i],x[1][0:i])
    #redDot.set_data(x[0][i], x[1][i])
    l.set_data(x[0][0:i],x[1][0:i])
    return l,
plt.title("Forced damped Oscillator")
plt.xlabel('t')
plt.ylabel('position')
myAnimation = animation.FuncAnimation(fig, animate, frames=np.arange(0, len(x[0]), 1), \
                                      interval=10, blit=True, repeat=True)
plt.show()


