import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as animation

data = np.genfromtxt("SolGeneral.cvs ",delimiter=",")

r1_x  = data[:,0]
r1_y  = data[:,1]
r1_z  = data[:,2]
r2_x  = data[:,3]
r2_y  = data[:,4]
r2_z  = data[:,5]
r3_x  = data[:,6]
r3_y  = data[:,7]
r3_z  = data[:,8]

fig = plt.figure()
ax = plt.axes(projection='3d')



def actualizar(i):
    ax.clear()
    ax.plot(r1_x [:i], r1_y[:i], r1_z[:i])
    ax.plot(r2_x [:i], r2_y[:i], r2_z[:i])
    ax.plot(r3_x [:i], r3_y[:i], r3_z[:i])


anim= animation.FuncAnimation(fig,actualizar,range(len(r1_x)), repeat = True)
anim.save('Tre.gif')
