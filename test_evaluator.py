from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import numpy as np

import matplotlib.pyplot as plt

# Read the data from the file
data = []
atoms = []
with open('evaluator.txt', 'r') as file:
    for line in file:
        values = line.strip().split()
        if len(values) == 3:
            atoms.append([float(values[0]), float(values[1]), float(values[2])])

        if len(values) == 4:
            data.append([float(values[0]), float(values[1]), float(values[2]), float(values[3])])

# Extract the coordinates and values
x = [row[0] for row in data]
y = [row[1] for row in data]
z = [row[2] for row in data]
alpha = [row[3] for row in data]
alpha = [val/max(alpha) for val in alpha]
size = [row[3]*80 for row in data]

# Extract the atom coordinates
atoms_x = [row[0] for row in atoms]
atoms_y = [row[1] for row in atoms]
atoms_z = [row[2] for row in atoms]

# Plot the data on a 3D graph
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x_scale=max(x)-min(x)
y_scale=max(y)-min(y)
z_scale=max(z)-min(z)

scale=np.diag([x_scale, y_scale, z_scale, 1.0])
scale=scale*(1.0/scale.max())
scale[3,3]=1.0

def short_proj():
  return np.dot(Axes3D.get_proj(ax), scale)

ax.get_proj=short_proj
ax.mouse_init()


ax.scatter(atoms_x, atoms_y, atoms_z, c='red', s=100)

# ax.scatter(x, y, z, alpha=alpha, c='black', s=size)
ax.scatter(x, y, z, c=alpha, cmap='viridis', s=20)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()