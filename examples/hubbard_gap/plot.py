import numpy as np
import matplotlib.pyplot as plt


m = np.genfromtxt("D_VS_MU.OUT").transpose()
x,y = m[0],m[1]
from scipy.interpolate import interp1d

from scipy.signal import savgol_filter

xs = np.linspace(np.min(x),np.max(x),len(x)*10)
ys = interp1d(x,y,fill_value=0.0,bounds_error=False)(xs)
#ys = savgol_filter(ys, 51, 5)
zs = np.gradient(ys)
zs2 = savgol_filter(zs, 51, 5)
plt.plot(xs,-zs)
plt.plot(xs,-zs2)
plt.show()
