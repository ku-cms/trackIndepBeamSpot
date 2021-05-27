# example_3d.py

import numpy as np
import matplotlib.pyplot as plt

def f(x, y):
    return np.sin(np.sqrt(x ** 2 + y ** 2))

def plot(plot_dir, plot_name):
    if plot_dir[-1] != "/":
        plot_dir += "/"
    
    x = np.linspace(-6, 6, 30)
    y = np.linspace(-6, 6, 30)
    
    X, Y = np.meshgrid(x, y)
    Z = f(X, Y)
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.contour3D(X, Y, Z, 50, cmap='binary')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    plt.savefig(plot_dir + plot_name + '.png', bbox_inches='tight')

plot("plots", "example_3d")

