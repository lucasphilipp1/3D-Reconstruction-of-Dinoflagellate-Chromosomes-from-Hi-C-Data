import alabtools, numpy as np

matrix=alabtools.analysis.get_simulated_hic("igm-model.hss")
array=matrix.toarray()
np.savetxt("simulated_hic.txt", array, fmt='%.5f')
