##############################################################################

import hicstraw
import matplotlib.pyplot as plt
import pandas as pd
import cooler
import cooltools
import numpy as np

#Smic1.0
#microadriaticum: 94 chromosomes

chrom_lengths_microadriaticum_Smic1_0=np.zeros(94+95)
chrom_names_microadriaticum_Smic1_0=[]

clr_microadriaticum_Smic1_0 = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files/GSE152150_Dino-HiC-cD1plus-R1.smic1.0.mapq_30.1000.mcool::resolutions/1000')

i=0
for chrom in clr_microadriaticum_Smic1_0.chromnames:
    chrom_lengths_microadriaticum_Smic1_0[i] = (clr_microadriaticum_Smic1_0.extent(chrom)[1]-clr_microadriaticum_Smic1_0.extent(chrom)[0])*1000
    chrom_names_microadriaticum_Smic1_0.append(clr_microadriaticum_Smic1_0.chromnames[i])
    i=i+1

chrom_length_ind_microadriaticum_Smic1_0 = np.argsort(chrom_lengths_microadriaticum_Smic1_0)
chrom_length_ind_microadriaticum_Smic1_0 = chrom_length_ind_microadriaticum_Smic1_0[::-1]

chrom_names_microadriaticum_Smic1_0 = [chrom_names_microadriaticum_Smic1_0[i] for i in chrom_length_ind_microadriaticum_Smic1_0]

chrom_lengths_microadriaticum_Smic1_0 = [chrom_lengths_microadriaticum_Smic1_0[i] for i in chrom_length_ind_microadriaticum_Smic1_0]

f, ax = plt.subplots(1,1)
plt.hist(chrom_lengths_microadriaticum_Smic1_0, bins = 100, range=(30000, max(chrom_lengths_microadriaticum_Smic1_0)))
plt.axvline(x=500000, color='orange')
plt.ylim((0, 25))
ax.set(
xlabel='scaffold length, bp',
ylabel='count')
plt.show()

#SMIC1.1N
#microadriaticum: 94 chromosomes

chrom_lengths_microadriaticum_Smic1_1N=np.zeros(94)
chrom_names_microadriaticum_Smic1_1N=[]

clr_microadriaticum_Smic1_1N = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files/GSE152150_Dino-HiC-cD1plus-R1.smic1.1N.mapq_30.1000.mcool::resolutions/1000')

i=0
for chrom in clr_microadriaticum_Smic1_1N.chromnames:
    chrom_lengths_microadriaticum_Smic1_1N[i] = (clr_microadriaticum_Smic1_1N.extent(chrom)[1]-clr_microadriaticum_Smic1_1N.extent(chrom)[0])*1000
    chrom_names_microadriaticum_Smic1_1N.append(clr_microadriaticum_Smic1_1N.chromnames[i])
    i=i+1

chrom_length_ind_microadriaticum_Smic1_1N = np.argsort(chrom_lengths_microadriaticum_Smic1_1N)
chrom_length_ind_microadriaticum_Smic1_1N = chrom_length_ind_microadriaticum_Smic1_1N[::-1]

chrom_names_microadriaticum_Smic1_1N = [chrom_names_microadriaticum_Smic1_1N[i] for i in chrom_length_ind_microadriaticum_Smic1_1N]

chrom_lengths_microadriaticum_Smic1_1N = [chrom_lengths_microadriaticum_Smic1_1N[i] for i in chrom_length_ind_microadriaticum_Smic1_1N]

f, ax = plt.subplots(1,1)
plt.hist(chrom_lengths_microadriaticum_Smic1_1N, bins = 100, range=(30000, max(chrom_lengths_microadriaticum_Smic1_1N)))
plt.axvline(x=500000, color='orange')
plt.ylim((0, 25))
ax.set(
xlabel='scaffold length, bp',
ylabel='count')
plt.show()

#minutum: 91 pseudochromosomes > 500kb
#longest being ∼11megabases (Mb) in size, with a median length of 6.7Mb

chrom_lengths_minutum = np.zeros(16003)
chrom_names_minutum = []

hic = hicstraw.HiCFile('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files/GSM4658994_L142+L533+L534.inter.hic')

i=0
for chrom in hic.getChromosomes():
  chrom_lengths_minutum[i] = chrom.length
  chrom_names_minutum.append(chrom.name)
  i=i+1

chrom_length_ind_minutum = np.argsort(chrom_lengths_minutum)
chrom_length_ind_minutum = chrom_length_ind_minutum[::-1]

chrom_names_minutum = [chrom_names_minutum[i] for i in chrom_length_ind_minutum]
chrom_lengths_minutum = [chrom_lengths_minutum[i] for i in chrom_length_ind_minutum]

#move ALL scaffold to end of all lists
chrom_length_ind_minutum = np.delete(chrom_length_ind_minutum, np.where(chrom_length_ind_minutum==chrom_names_minutum.index('ALL')))
chrom_length_ind_minutum = np.append(chrom_length_ind_minutum, chrom_names_minutum.index('ALL'))

chrom_names_minutum = [chrom_names_minutum[i] for i in chrom_length_ind_minutum]
chrom_lengths_minutum = [chrom_lengths_minutum[i] for i in chrom_length_ind_minutum]

f, ax = plt.subplots(1,1)
plt.hist(chrom_lengths_minutum, bins = 100, range=(30000, max(chrom_lengths_minutum)))
plt.axvline(x=500000, color='orange')
plt.ylim((0, 25))
ax.set(
xlabel='scaffold length, bp',
ylabel='count')
plt.show()

#print(len([i for i in chrom_lengths_minutum if i > 500000]))
#print([chrom_names_minutum[chrom_lengths_minutum.index(i)] for i in chrom_lengths_minutum if i > 500000])
#print(chrom_names_minutum[91])

##############################################################################