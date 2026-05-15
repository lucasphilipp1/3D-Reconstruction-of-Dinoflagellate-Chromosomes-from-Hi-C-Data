import alabtools
import numpy as np

coords = alabtools.HssFile("igm-model.hss", 'r')
array = coords.get_coordinates()

for i in range(array.shape[1]):
    filename = f"coordinates_{i}.txt"
    np.savetxt(filename, array[:, i, :], fmt="%.6f", delimiter="\t")
