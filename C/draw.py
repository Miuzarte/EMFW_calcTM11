import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('Ez.csv', header=None)
matrix = df.values
x, y = np.mgrid[:matrix.shape[0], :matrix.shape[1]]

fig = plt.figure(figsize=(12, 6))

ax1 = fig.add_subplot(121, projection='3d')
surf = ax1.plot_surface(x, y, matrix, cmap='viridis')
fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=5)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('TM11 Ez Surface')

ax2 = fig.add_subplot(122)
contour = ax2.contourf(x, y, matrix, cmap='viridis')
fig.colorbar(contour, ax=ax2, shrink=0.5, aspect=5)
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_title('TM11 Ez Contour')

plt.show()
