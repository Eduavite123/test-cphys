import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt('ham_prima.dat')
x=data[:,0]
y=data[:,1]

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(x,y)
#ax.set_ylim(ymin=-1.3e10,ymax=-1.5e10)

plt.xlabel('Tiempo (horas)')
plt.ylabel('$H^\prime=H-\omega p_{\phi}$')
fig.set_tight_layout(True)

plt.show()

file_out='ham_prima'
fig.savefig("{}.png".format(file_out))
