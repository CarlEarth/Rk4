import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
TWOPI = 2*np.pi

fig, ax = plt.subplots()
#xy1 = np.array([0.2, 0.2])
xy2 = np.array([1.2, -0.1])
#circle = patches.Circle(xy1, 0.1, color = "r")
rec=patches.Rectangle(xy2,0.2,0.2, color = "r")
ax.add_patch(rec)
ax = plt.axis([-1,1,-1,1])

Y = np.loadtxt('time_position')
x=Y.T
"""
fig, ax = plt.subplots()
ax = plt.axis([-1,1,-1,1])
#l, =plt.plot(x[0][0:1],x[1][0:1])
xy2 = np.array([0.5, 0.5])
#rec =patches.Rectangle(xy2,1,1, color = "r")
circle = patches.Circle(xy2, 0.1, color = "r")
ax.add_patch(circle)
#def animate(i):
    #rec.set_x(x[1][i]-0.5)
#return rec,
"""
def animate(i):
    rec.set_x(x[1][i]-0.1)
    return rec,
plt.title("Forced damped Oscillator")
plt.xlabel('Position (m)')
plt.ylabel('')
myAnimation = animation.FuncAnimation(fig, animate, frames=np.arange(0, len(x[0]), 1), \
                                      interval=10, blit=True, repeat=True)
plt.show()


