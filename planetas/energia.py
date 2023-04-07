import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt('energia.dat')
x=data[:,0]
y=data[:,1]

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(x,y)
ax.set_ylim(ymin=-0.003,ymax=-0.002)

plt.xlabel('Tiempo')
plt.ylabel('Energía')
fig.set_tight_layout(True)

plt.show()

file_out='energia'
fig.savefig("{}.pdf".format(file_out))
