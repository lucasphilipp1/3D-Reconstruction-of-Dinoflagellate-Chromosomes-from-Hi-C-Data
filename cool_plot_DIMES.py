import os
import cooler
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import scipy
import scipy.spatial
from matplotlib.ticker import EngFormatter
from matplotlib.colors import LogNorm

bp_formatter = EngFormatter('b')

#Allows us to plot the data points on a log scale
norm = LogNorm(vmax=500)


def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

cur_path = os.path.dirname(__file__)


#Example cool file that works with the plot 
file_path_good = r"C:\Users\alyss\Desktop\Dinoflagellate\HIPPS-DIMES\data\Rao2014-GM12878-MboI-allreps-filtered.1000kb.cool"

#DIMES control file 
file_path = r"C:\Users\alyss\Desktop\Dinoflagellate\HIPPS-DIMES\data\GSM2745898_17MAY11_HSC02RPAB-PE50_CCHiC-HeLa-NS-R1_23-05-2011_kittlere.5.10000.cool"

clr = cooler.Cooler(file_path)

chromstarts = []
for i in clr.chromnames:
    chromstarts.append(clr.extent(i)[0])

f, ax = plt.subplots(figsize=(7,6))

#Specify the chromosome number here (applies to rest of the code)
chromosome_num = "17"

im = ax.matshow(
    clr.matrix(balance=False).fetch('chr'+chromosome_num),
    norm=norm,
    extent=(0,clr.chromsizes['chr'+chromosome_num], clr.chromsizes['chr'+chromosome_num], 0)
);
plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='raw counts');
ax.set_title('Chromosome ' + chromosome_num, y=1.08)
ax.set_ylabel('position, Mb')
format_ticks(ax)


