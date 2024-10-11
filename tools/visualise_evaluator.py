from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import numpy as np

import matplotlib.pyplot as plt

# Read the data from the file
data = []
atoms = []
with open('viability.dat', 'r') as file:
    for line in file:
        values = line.strip().split()
        if len(values) == 3:
            atoms.append([float(values[0]), float(values[1]), float(values[2])])

        if len(values) == 4:
            data.append([float(values[0]), float(values[1]), float(values[2]), float(values[3])])

# Drop all rows with row[3] less than 0.2
# data = [row for row in data if row[3] > 0.4] # for 4-body only
max_val = max([row[3] for row in data])
avg_val = sum([row[3] for row in data])/len(data)
print(f"Max value: {max_val}")
print(f"Avg value: {avg_val}")

# Extract the coordinates and values
x = [row[0] for row in data]
y = [row[1] for row in data]
z = [row[2] for row in data]
# alpha = [row[3] for row in data]
alpha = [(row[3]*100)*2 for row in data]
alpha = [val/max(alpha) for val in alpha]
# size = [(row[3]*10)**3 for row in data] # for 2-body only
# size = [(row[3]*100)*2 for row in data] # for 2-body and 3-body
# size = [row[3] for row in data] # for 4-body only
# size = [(row[3]*2000) for row in data] # for all
size = [(row[3]*1000) for row in data] # for all

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


ax.scatter(atoms_x, atoms_y, atoms_z, c='red', alpha=1.0, s=200)

# ax.scatter(x, y, z, alpha=alpha, c='black', s=size)
ax.scatter(x, y, z, c=alpha, cmap='viridis', s=size)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
