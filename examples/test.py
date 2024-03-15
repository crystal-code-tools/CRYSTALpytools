import matplotlib.pyplot as plt

import CRYSTALpytools.plot as cfplt
from CRYSTALpytools.crystal_io import Properties_output

data = Properties_output().read_electron_dos('data/doss_96.DOSS')
print(type(data.energy[0]))
cfplt.plot_electron_dos(data)

plt.show()
