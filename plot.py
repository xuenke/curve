import csv
import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt(open('fit.csv', 'rb'), delimiter=",", skiprows = 0)
b = np.loadtxt(open('fit2.csv', 'rb'), delimiter=",", skiprows = 0)

plt.plot(a[:,0], a[:, 1], 'ob-')
plt.plot(b[:,0], b[:, 1], '.r-')
plt.axis("equal")
plt.show()