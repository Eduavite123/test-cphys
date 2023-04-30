import matplotlib.pyplot as plt
import numpy as np
import io

data=np.loadtxt("norma.txt")

fig=plt.plot(data, color='red', label='Data')
plt.xlabel('Tiempo n')
plt.ylabel('$||\psi||^2$')
plt.title('$\lambda$=0.1')

plt.show()