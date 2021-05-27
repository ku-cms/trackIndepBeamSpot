# example_3d.py

# Matplotlib 3D plot examples:
# https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html

import numpy as np
import matplotlib.pyplot as plt

def f(x, y):
    return np.sin(np.sqrt(x ** 2 + y ** 2))

def plot(plot_dir, plot_name):
    if plot_dir[-1] != "/":
        plot_dir += "/"
    
    x = np.linspace(-6, 6, 30)
    y = np.linspace(-6, 6, 30)
    
    # X, Y, and Z are 2D matrices
    X, Y = np.meshgrid(x, y)
    Z = f(X, Y)

    #print("X = {0}".format(X))
    #print("Y = {0}".format(Y))
    #print("Z = {0}".format(Z))

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    title = "contour"
    file_name = "{0}{1}_{2}.png".format(plot_dir, plot_name, title)
    ax.contour3D(X, Y, Z, 50, cmap='binary')
    ax.set_title(title)
    plt.savefig(file_name, bbox_inches='tight')

    title = "wireframe"
    file_name = "{0}{1}_{2}.png".format(plot_dir, plot_name, title)
    ax.plot_wireframe(X, Y, Z, color='black')
    ax.set_title(title)
    plt.savefig(file_name, bbox_inches='tight')
    
    title = "surface"
    file_name = "{0}{1}_{2}.png".format(plot_dir, plot_name, title)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.set_title(title)
    plt.savefig(file_name, bbox_inches='tight')

plot("plots", "example_3d")

