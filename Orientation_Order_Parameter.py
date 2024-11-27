#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:58:52 2023

code adapted from:
https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.248101

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

def compute_Orientation_OP(xyz, chrom_start, chrom_end, vec_length):
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

def compute_FFT_from_Oij(Oijy, resolution):
    from scipy.fftpack import fft, fftfreq
 
    N = np.shape(Oijy)[0]
    xf = fftfreq(N, resolution)[:N//2]
    yf = fft(Oijy)
    return xf, 2.0/N * np.abs(yf[0:N//2])

def tolerant_mean(arrs):
    lens = [len(i) for i in arrs]
    arr = np.ma.empty((np.max(lens),len(arrs)))
    arr.mask = True
    for idx, l in enumerate(arrs):
        arr[:len(l),idx] = l
    return arr.mean(axis = -1), arr.std(axis=-1)

Oijx_list=[]
Oijy_list=[]

xf_list=[]
yf_list=[]

res=5000 #Hi-C resolution
cutoff=5*10**6 #ignores high frequency noise at large primary sequence separation, and allows for differently sized chromosomes to be superimposed

Oijx_list.clear()
Oijy_list.clear()

xf_list.clear()
yf_list.clear()

#calculate orientation order parameter (OOP) and its Fourier transform
for i in range(1, 76):
    #df = np.loadtxt('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/CSynth 3D smicroadriaticum/structures/combined/symbiodinium_microadriaticum_chr'+str(i)+'_3D.xyz')
    df = np.loadtxt('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/CSynth 3D skawagutii/structures/s_kawagutii_V3_HiC_scaffold_'+str(i)+'.xyz')
    #df = np.loadtxt('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/GSE18199_Globules/equilibrium/equilibrium'+str(i)+'.dat', skiprows=1)    
        
    structure = df[:,1:] #remove first column, don't confuse primary sequence with spatial position
    #structure = df
    Oijx, Oijy = compute_Orientation_OP(xyz=structure, chrom_start=0, chrom_end=df.shape[0], vec_length=4)
    
    Oijx=Oijx*res
    
    Oijy=np.delete(Oijy,range(-round((max(Oijx)-cutoff)/res),0)) #ignore high-frequency data at large primary sequence separation
    Oijx=np.delete(Oijx,range(-round((max(Oijx)-cutoff)/res),0)) #ignore high-frequency data at large primary sequence separation
    
    Oijx_list.append(Oijx)
    Oijy_list.append(Oijy)
    
    xf, yf = compute_FFT_from_Oij(Oijy=Oijy, resolution=res)
    xf_list.append(xf)
    yf_list.append(yf)
    print(i)


#plot OOP
for i in range(len(yf_list)):
    #plt.plot(Oijx_list[i],Oijy_list[i], color='#1f77b4', alpha = 0.04) #microadriaticum
    plt.plot(Oijx_list[i],Oijy_list[i], color='#2ca02c', alpha = 0.04) #kawagutii
    #plt.plot(Oijx_list[i],Oijy_list[i], color='#bcbd22', alpha = 0.04) #equilibrium globule
    #plt.plot(Oijx_list[i],Oijy_list[i], color='#000000', alpha = 1) #CLC
    
plt.xscale('log',base=10)
plt.xlabel('bp', fontsize=16)
#plt.ylabel('Orientation Order Parameter', fontsize=16)
#plt.title('Symbiodinium microadriaticum')
#plt.title('Symbiodinium kawagutii')
#plt.title('Equilibrium Globlues')
ax = plt.gca()
ax.set_ylim([-1, 1])
ax.tick_params(axis='both', which='major', labelsize=16)
#ax.errorbar(Oijx_avg, Oijy_avg, yerr=error[0:1001])

#plt.axvline(x=7.3418*10**4, color='#2ca02c', linestyle='dashed') #mean tandem gene array length (kawagutii)
#ax.axvspan(5.4902*10**4, 9.1933*10**4, alpha=0.5, color='#2ca02c') #+/- one std of tandem gene array length (kawagutii)

plt.axvline(x=5.8961*10**4, color='#1f77b4', linestyle='dashed') #mean tandem gene array length (microadriaticum)
ax.axvspan(4.4072*10**4, 7.3851*10**4, alpha=0.5, color='#1f77b4') #+/- one std of tandem gene array length (microadriaticum)

#plt.axvline(x=2.3123*10**6 , color='#2ca02c', linestyle='dotted') #mean tad size (kawagutii)
plt.axvline(x=1.6559*10**6 , color='#1f77b4', linestyle='dotted') #mean tad size (microadriaticum)

#plt.axvline(x=4.3748*10**3, color='#2ca02c', linestyle='dashdot') #mean gene length (kawagutii)
#plt.axvline(x=5.0781*10**3 , color='#1f77b4', linestyle='dashdot') #mean gene length (microadriaticum)

plt.show()

#plot OOP Fourier transform
for i in range(len(yf_list)):
    #plt.plot(xf_list[i],yf_list[i], color='#1f77b4', alpha = 0.04) #microadriaticum
    plt.plot(xf_list[i],yf_list[i], color='#2ca02c', alpha = 0.04) #kawagutii
    #plt.plot(xf_list[i],yf_list[i], color='#bcbd22', alpha = 0.04) #equilibrium globule
    #plt.plot(xf_list[i],yf_list[i], color='#000000', alpha = 1) #CLC
    
#plt.plot(xf,yf)
plt.xlabel('1/bp', fontsize=16)
#plt.ylabel('FT of Orientation Order Parameter', fontsize=12)
plt.xscale('log',base=10) 

plt.axvline(x=1/(7.3418*10**4), color='#2ca02c', linestyle='dashed') #mean tandem gene array length (kawagutii)
plt.axvline(x=1/(1.5122*10**6) , color='#2ca02c', linestyle='dotted') #mean tad size (kawagutii)
#plt.axvline(x=1/(4.3748*10**3), color='#2ca02c', linestyle='dashdot') #mean gene length (kawagutii)

#plt.axvline(x=1/(5.8961*10**4), color='#1f77b4', linestyle='dashed') #mean tandem gene array length (microadriaticum)
#plt.axvline(x=1/(2.7436*10**6), color='#1f77b4', linestyle='dotted') #mean tad size (microadriaticum)
#plt.axvline(x=1/(5.0781*10**3), color='#1f77b4', linestyle='dashdot') #mean gene length (microadriaticum)

#plt.axvline(x=1/(2**10**5), color='#000000', linestyle='dashed')
#plt.axvline(x=1/(2*0.75*10**6), color='#000000', linestyle='dashed')
#plt.title('symbiodinium_microadriaticum')
#plt.title('symbiodinium_kawagutii')
#plt.title('Equilibrium Globlues')
ax = plt.gca()
ax.set_ylim([0, 0.05])
ax.tick_params(axis='both', which='major', labelsize=16)

#ax.axvspan(1/(7.3851*10**4), 1/(4.4072*10**4), alpha=0.5, color='#1f77b4') #+/- one std of tandem gene array length (microadriaticum)
ax.axvspan(1/(9.1933*10**4), 1/(5.4902*10**4), alpha=0.5, color='#2ca02c') #+/- one std of tandem gene array length (kawagutii)

plt.show()



   

