# example_2d.py

# Plot a Basic 2D Histogram using Matplotlib
# https://www.python-graph-gallery.com/83-basic-2d-histograms-with-matplotlib

import matplotlib.pyplot as plt
import numpy as np

# create data
#x = np.random.normal(size=50000)
#y = 2 * x + np.random.normal(size=50000)

x = [1.0, 2.0]
y = [1.0, 2.0]

# Big bins
#plt.hist2d(x, y, bins=(100, 100), cmap=plt.cm.jet)
plt.hist2d(x, y, bins=(10, 10), cmap=plt.cm.viridis)
plt.show()

