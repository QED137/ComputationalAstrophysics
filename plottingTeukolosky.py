
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("frame_200.csv")
Nr = 250
Nth = 180

psi = df['re'].values.reshape(Nr, Nth)
plt.imshow(psi, aspect='auto', cmap='inferno')
plt.colorbar()
plt.show()
