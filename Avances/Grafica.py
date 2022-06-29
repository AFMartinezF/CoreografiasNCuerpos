import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


data = np.genfromtxt('Resultados.csv', delimiter=',')
r1_x  = data[1:,0]
r1_y  = data[1:,1]
r2_x  = data[1:,2]
r2_y  = data[1:,3]
r3_x  = data[1:,4]
r3_y  = data[1:,5]

fig = plt.figure()
ax = fig.gca()

def actualizar(i):
    ax.clear()
    ax.plot(r1_x [:i], r1_y[:i])
    ax.plot(r2_x [:i], r2_y[:i])
    ax.plot(r3_x [:i], r3_y[:i])

anim= animation.FuncAnimation(fig,actualizar,range(len(r1_x)), repeat = True)
anim.save('Grafica.gif')
