from matplotlib import pyplot as plt

import numpy as np

data = np.loadtxt('basis_functions.txt')
for i in range(1, data.shape[1]):
    plt.plot(data[:, 0], data[:, i], label=f'N_{i}')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('plot.pdf')
