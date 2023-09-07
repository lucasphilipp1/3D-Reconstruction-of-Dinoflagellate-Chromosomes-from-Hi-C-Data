#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:58:52 2023

@author: lucasphilipp
"""

import numpy as np
import matplotlib.pyplot as plt

# Copyright (c) 2020-2023 The Center for Theoretical Biological Physics (CTBP) - Rice University
# This file is from the Open-MiChroM project, released under the MIT License.

R"""
The :class:`~.cndbTools` class perform analysis from **cndb** or **ndb** - (Nucleome Data Bank) file format for storing an ensemble of chromosomal 3D structures.
Details about the NDB/CNDB file format can be found at the `Nucleome Data Bank <https://ndb.rice.edu/ndb-format>`__.
"""

#########################################################################################
#### Analysis start here!
#########################################################################################

def compute_Orientation_OP(xyz, chrom_start=0, chrom_end=1000, vec_length=4):
    from collections import OrderedDict
    import itertools
    R"""
    Calculates the Orientation Order Parameter OP. Details are decribed in "Zhang, Bin, and Peter G. Wolynes. "Topology, structures, and energy landscapes of human chromosomes." Proceedings of the National Academy of Sciences 112.19 (2015): 6062-6067."

    Args:
        xyz (:math:`(frames, beadSelection, XYZ)` :class:`numpy.ndarray`, required):
            Array of the 3D position of the selected beads for different frames extracted by using the :code: `xyz()` function.
        chrom_start (int, required):
            First bead to consider in the calculations (Default value = 0).
        chrom_end (int, required):
            Last bead to consider in the calculations (Default value = 1000).
        vec_length (int, required):
            Number of neighbor beads to build the vector separation :math:`i` and :math:`i+4` if vec_length is set to 4. (Default value = 4).


    Returns:
        Oijx:class:`numpy.ndarray`:
            Returns the genomic separation employed in the calculations.
        Oijy:class:`numpy.ndarray`:
            Returns the Orientation Order Parameter OP as a function of the genomic separation.
    """

    vec_rij = []
    for i in range(chrom_start, chrom_end-vec_length):
        # from a trajectory, gets the vector ri,i+vec_length
        vec_rij.append(xyz[i+vec_length]-xyz[i])

    dot_ri_rj = []
    ij = []
    for i in itertools.combinations_with_replacement(range(1, chrom_end-vec_length), 2):
        # dot product between all vector r,r+_vec_length
        dot_ri_rj.append(np.dot(
            vec_rij[i[0]]/np.linalg.norm(vec_rij[i[0]]), vec_rij[i[1]]/np.linalg.norm(vec_rij[i[1]])))
        ij.append(i[1]-i[0])  # genomic separation

    d = OrderedDict()
    for k, v in zip(ij, dot_ri_rj):
        # sum the values v for each genomic separation k
        d[k] = d.get(k, 0) + v

    Oijx = []
    Oijy = []
    for i in range(1, np.size(list(d.keys()))+1):
        Oijx.append(list(d.keys())[i-1]+1)
        # gets Oij normalized by the number of elements. For example, in a chromosome of length 100, genomic distance 1 has more elements (99) considered than genomic distance 100 (1).
        Oijy.append(list(d.values())[i-1]/(chrom_end-list(d.keys())[i-1]))

    return np.asarray(Oijx), np.asarray(Oijy)

def compute_FFT_from_Oij(Oijy, lowcut=1, highcut=500, order=5):
    from scipy.fftpack import fft
    R"""
    Calculates the Fourier transform of the Orientation Order Parameter OP. Details are decribed in "Zhang, Bin, and Peter G. Wolynes. "Topology, structures, and energy landscapes of human chromosomes." Proceedings of the National Academy of Sciences 112.19 (2015): 6062-6067."

    Args:
        xyz (:math:`(frames, beadSelection, XYZ)` :class:`numpy.ndarray`, required):
            Array of the 3D position of the selected beads for different frames extracted by using the :code: `xyz()` function.
        lowcut (int, required):
            Filter to cut low frequencies (Default value = 1).
        highcut (int, required):
            Filter to cut high frequencies (Default value = 500).
        order (int, required):
            Order of the Butterworth filter obtained from `scipy.signal.butter <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html>`__. (Default value = 5).


    Returns:
        xf:class:`numpy.ndarray`:
            Return frequencies.
        yf:class:`numpy.ndarray`:
            Returns the Fourier transform of the Orientation Order Parameter OP in space of 1/Chrom_Length.
    """
    N = np.shape(Oijy)[0]
    xf = np.linspace(1, N//2, N//2)
    yf = fft(Oijy)/len(Oijy)
    return (xf[0:N//2]-1)/N, np.abs(yf[0:N//2])


def tolerant_mean(arrs):
    lens = [len(i) for i in arrs]
    arr = np.ma.empty((np.max(lens),len(arrs)))
    arr.mask = True
    for idx, l in enumerate(arrs):
        arr[:len(l),idx] = l
    return arr.mean(axis = -1), arr.std(axis=-1)

Oijx_list=[]
Oijy_list=[]

for i in range(1, 101):
#for i in range(1, 2):
    #df = np.loadtxt('/Users/lucasphilipp/Downloads/bminutum_pseudochromosome_'+str(i)+'_CF_60_SP_0_PP-4_3D.xyz')
    df = np.loadtxt('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_3D_structures/CSynth 3D bminutum/structures/bminutum_pseudochromosome_'+str(i)+'_3D.xyz')
    structure = df[:,1:] #remove first column
    Oijx, Oijy = compute_Orientation_OP(xyz=structure, chrom_start=0, chrom_end=df.shape[0], vec_length=4)
    
    Oijx_list.append(Oijx)
    Oijy_list.append(Oijy)
    print(i)

curr_max = np.zeros(1)
#get longest Oijx
for o in Oijx_list:
    if o.shape[0]>curr_max.shape[0]:
        curr_max=o

#average Oijy at each frequency
Oijy_avg, error = tolerant_mean(Oijy_list)

Oijx_avg=np.arange(0,curr_max.shape[0])
Oijx_avg=Oijx_avg*5000

Oijy_avg=np.delete(Oijy_avg,range(-round((max(Oijx_avg)-5*10**6)/5000),0)) #ignore high-frequency data at large separation
Oijx_avg=np.delete(Oijx_avg,range(-round((max(Oijx_avg)-5*10**6)/5000),0)) #ignore high-frequency data at large separation

plt.plot(Oijx_avg,Oijy_avg)
plt.xscale('log',base=10)
plt.xlabel('bp')
plt.ylabel("Orientation Order Parameter")
plt.title('breviolum_minutum')
plt.show()

xf, yf = compute_FFT_from_Oij(Oijy=Oijy_avg)
plt.plot(xf/5000,yf)
plt.xscale('log',base=10) 
plt.xlabel('1/bp')
plt.ylabel('FT of Orientation Order Parameter')
plt.title('breviolum_minutum')
plt.show()

Oijx_list.clear()
Oijy_list.clear()

#for i in range(1, 95):
for i in range(1, 2):
    #df = np.loadtxt('/Users/lucasphilipp/Downloads/symbiodinium_microadriaticum_coccoid_chr'+str(i+2)+'_CF_60_SP_0_PP-4_3D.xyz')
    df = np.loadtxt('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_3D_structures/CSynth 3D smicroadriaticum/structures/symbiodinium_microadriaticum_coccoid_chr'+str(i)+'_3D.xyz')
    structure = df[:,1:] #remove first column
    Oijx, Oijy = compute_Orientation_OP(xyz=structure, chrom_start=0, chrom_end=df.shape[0], vec_length=4)
    
    Oijx_list.append(Oijx)
    Oijy_list.append(Oijy)
    print(i)

curr_max = np.zeros(1)
#get longest oijx
for o in Oijx_list:
    if o.shape[0]>curr_max.shape[0]:
        curr_max=o

#average Oijy at each frequency
Oijy_avg, error = tolerant_mean(Oijy_list)

Oijx_avg=np.arange(0,curr_max.shape[0])
Oijx_avg=Oijx_avg*5000

testy=Oijy_avg
testx=Oijx_avg

Oijy_avg=testy
Oijx_avg=testx

Oijy_avg=np.delete(Oijy_avg,range(-round((max(Oijx_avg)-2*10**7)/5000),0)) #ignore high-frequency data at large separation
Oijx_avg=np.delete(Oijx_avg,range(-round((max(Oijx_avg)-2*10**7)/5000),0)) #ignore high-frequency data at large separation

plt.plot(Oijx_avg,Oijy_avg)
plt.xscale('log',base=10)
plt.xlabel('bp')
plt.ylabel("Orientation Order Parameter")
plt.title('symbiodinium_microadriaticum_coccoid')
plt.show()

xf, yf = compute_FFT_from_Oij(Oijy=Oijy_avg)
plt.plot(xf/5000,yf)
plt.xscale('log',base=10) 
plt.xlabel('1/bp')
plt.ylabel('FT of Orientation Order Parameter')
plt.title('symbiodinium_microadriaticum_coccoid')
plt.show()

#/Users/lucasphilipp/Downloads/cholesteric_CSynth_D4.txt


   

