import numpy as np

m = np.genfromtxt(".mpsfolder/MEASURE_SZ.OUT").transpose()

print("S = ",np.sum(m[1]))
