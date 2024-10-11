# %%
import numpy as np
from scipy.interpolate import griddata
from mayavi import mlab
# from mayavi.api import Engine



# %%
a = 1.0
b = 1.0
c = 1.0
n_x = 0
n_y = 0
n_z = 0

# Read the data from the file
data = []
atoms = []
with open('../viability.dat', 'r') as file:
    for line in file:
        # if line starts with #grid, read the three values as n_x, n_y, n_z
        if line.startswith("#grid"):
            values = line.strip().split()
            n_x = int(values[1])
            n_y = int(values[2])
            n_z = int(values[3])
            continue
        # if line starts with #lat, read the three values as a, b, c
        if line.startswith("#lat"):
            values = line.strip().split()
            a = float(values[1])
            b = float(values[2])
            c = float(values[3])
            continue
        values = line.strip().split()
        if len(values) == 3:
            atoms.append([float(values[0]), float(values[1]), float(values[2])])
        if len(values) == 4:
            data.append([float(values[0]), float(values[1]), float(values[2]), float(values[3])])


# Extract the coordinates and values
x = np.array([row[0] for row in data])
y = np.array([row[1] for row in data])
z = np.array([row[2] for row in data])
values = np.array([row[3] for row in data])

# # scale the data positions
# x = x * a
# y = y * b
# z = z * c

# scale atom locations
atoms = np.array(atoms)
atoms[:,0] = atoms[:,0] * n_x
atoms[:,1] = atoms[:,1] * n_y
atoms[:,2] = atoms[:,2] * n_z


# %%
# Create a 3D grid
grid_x, grid_y, grid_z = np.mgrid[x.min():x.max():complex(n_x), y.min():y.max():complex(n_y), z.min():z.max():complex(n_z)]


# %%
# Interpolate data onto the 3D grid
grid_values = griddata((x, y, z), values, (grid_x, grid_y, grid_z), method='nearest')
grid_values = np.nan_to_num(grid_values)

# %%
# Function to plot atoms as spheres
def plot_atoms_as_spheres(ax, atoms, radius=0.5, resolution=(10, 5)):
    u = np.linspace(0, 2 * np.pi, resolution[0])
    v = np.linspace(0, np.pi, resolution[1])
    for atom in atoms:
        atom_x, atom_y, atom_z = atom
        # Generate sphere coordinates
        x_sphere = atom_x + radius * np.outer(np.cos(u), np.sin(v))
        y_sphere = atom_y + radius * np.outer(np.sin(u), np.sin(v))
        z_sphere = atom_z + radius * np.outer(np.ones(np.size(u)), np.cos(v))
        # Plot spheres with lower resolution
        ax.plot_surface(x_sphere, y_sphere, z_sphere, color='red', alpha=1.0, rstride=1, cstride=1)
        
def plot_atoms_as_scatter(ax, atoms, size=100, color='red'):
    ax.scatter(atoms[:,0], atoms[:,1], atoms[:,2], s=size, c=color, alpha=1.0)

# Get atom bonds within radius
def get_bonds(atoms, grid, abc, radius=1.0):
    bonds = []
    for i, atom1 in enumerate(atoms):
        for j, atom2 in enumerate(atoms):
            if i != j:
                distance = np.linalg.norm( ( atom1 - atom2 ) * abc / grid)
                if distance <= radius:
                    bonds.append((i, j))
    return bonds

# Plot the atoms as spheres
def plot_atoms(atoms, bonds, radius=0.1, color=(1, 0, 0)):
    x_list = []
    y_list = []
    z_list = []
    for atom in atoms:
        x_sphere, y_sphere, z_sphere = atom
        x_list.append(x_sphere)
        y_list.append(y_sphere)
        z_list.append(z_sphere)
        # mlab.points3d(x_sphere, y_sphere, z_sphere, scale_factor=radius, color=color)
    pts = mlab.points3d(x_list, y_list, z_list, scale_factor=radius, color=color)

    pts.mlab_source.dataset.lines = np.array(bonds)

    tube = mlab.pipeline.tube(pts, tube_radius=1.0)
    tube.filter.radius_factor = 1.
    mlab.pipeline.surface(tube, color=color)


# # Draw the lattice bounding box
# def draw_bounding_box():
#     mlab.plot3d([0, n_x, n_x, 0, 0, 0, n_x, n_x, 0, 0, 0, 0, 0, n_x],
#                  [0, 0, n_y, n_y, 0, 0, 0, 0, 0, n_y, n_y, 0, 0, 0],
#                  [0, 0, 0, 0, 0, n_z, n_z, n_z, n_z, n_z, n_z, n_z, 0, 0],
#                  tube_radius=0.01, color=(0, 0, 0))  # Change the color as needed
    
    # Draw the lattice bounding box
def draw_bounding_box(a, b, c):
    # Define the vertices of the bounding box
    vertices = np.array([[0.0, 0.0, 0.0], [a, 0.0, 0.1], [a, b, 0.0], [0.0, b, 0.0],
                         [0.0, 0.0, c], [a, 0.0, c], [a, b, c], [0.0, b, c]])
    
    # Draw the lines connecting the vertices to create the box
    edges = [
        [vertices[0], vertices[1]], [vertices[1], vertices[2]], [vertices[2], vertices[3]], [vertices[3], vertices[0]],  # Bottom face
        [vertices[4], vertices[5]], [vertices[5], vertices[6]], [vertices[6], vertices[7]], [vertices[7], vertices[4]],  # Top face
        [vertices[0], vertices[4]], [vertices[1], vertices[5]], [vertices[2], vertices[6]], [vertices[3], vertices[7]],  # Vertical edges
    ]
    
    connections = ((0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 5), (2, 6), (3, 7), (4, 5), (5, 6), (6, 7), (7, 4), (4, 0), (5, 1), (6, 2), (7, 3))
    x_list = []
    y_list = []
    z_list = []
    for vertex in vertices:
        x_list.append(vertex[0])
        y_list.append(vertex[1])
        z_list.append(vertex[2])
    pts = mlab.points3d(x_list, y_list, z_list, scale_factor=0.1, color=(0, 0, 0))
    pts.mlab_source.dataset.lines = np.array(connections)
    tube = mlab.pipeline.tube(pts, tube_radius=0.5)
    tube.filter.radius_factor = 1.
    mlab.pipeline.surface(tube, color=(0, 0, 0))

    # for edge in edges:
    #     mlab.plot3d(*zip(*edge), tube_radius=0.01, color=(0, 0, 0))


# %%
print(f"X: {x.min()} - {x.max()}")
print(f"Y: {y.min()} - {y.max()}")
print(f"Z: {z.min()} - {z.max()}")
print(f"gX: {grid_x.min()} - {grid_x.max()}")


cont3d = mlab.contour3d(grid_values, contours=10, transparent=True)
# ax1 = mlab.axes( color=(1,1,1), nb_labels=4 )

# Call the function to plot atoms
bonds = get_bonds(atoms, [n_x, n_y, n_z], [a, b, c], radius=2.5)
plot_atoms(atoms, bonds, radius=10.0)

# Call the function to draw the bounding box
draw_bounding_box(n_x, n_y, n_z)
# draw_bounding_box()

mlab.show()
