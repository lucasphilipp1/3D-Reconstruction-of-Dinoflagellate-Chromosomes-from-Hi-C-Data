"""
a minor modification of:
https://cooltools.readthedocs.io/en/latest/notebooks/contacts_vs_distance.html
"""

import cooler
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
from matplotlib.colors import LogNorm
from cooltools.lib.numutils import adaptive_coarsegrain, interp_nan

bp_formatter = EngFormatter('b') #label in bp not kilobase-pairs or megabase-pairs

#logarithmic color scale
norm = LogNorm(vmax=0.02)


def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)


clr = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/GSE152150_HiC-Dplus.smic1.1N.mapq_30.1000.mcool::/resolutions/5000')

chromstarts = []
for i in clr.chromnames:
    chromstarts.append(clr.extent(i)[0])

f, ax = plt.subplots(figsize=(7,6))

#specific the name of the chromosome here:
#microadriaticum: chr#_pilon
#kawagutii: HiC_scaffold_#
chromosome_num = "chr1_pilon"

im = ax.matshow(
    clr.matrix(balance=True).fetch(chromosome_num),
    norm=norm,
    extent=(0,clr.chromsizes[chromosome_num], clr.chromsizes[chromosome_num], 0),
    
);
plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='raw counts');
ax.set_title('Chromosome ' + chromosome_num, y=1.08)
ax.set_ylabel('position, Mb')
format_ticks(ax)

region = (chromosome_num, 0, clr.chromsizes[chromosome_num])
extent=(0,clr.chromsizes[chromosome_num], clr.chromsizes[chromosome_num], 0),

cg = adaptive_coarsegrain(clr.matrix(balance=True).fetch(region),
                              clr.matrix(balance=False).fetch(region),
                              cutoff=3, max_levels=8)

cgi = interp_nan(cg)

f, axs = plt.subplots(
    figsize=(18,5),
    nrows=1,
    ncols=1,
    sharex=True, sharey=True)

im3 = axs.matshow(cgi, norm=norm, extent=extent)

#format_ticks(axs, rotate=False)
plt.gca().xaxis.tick_bottom()

plt.colorbar(im3, ax=axs, fraction=0.046, label='Normalized Contact Frequency')

clr = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/Symbiodinium_kawagutii_5000.cool')

chromstarts = []
for i in clr.chromnames:
    chromstarts.append(clr.extent(i)[0])

f, ax = plt.subplots(figsize=(7,6))

#specific the name of the chromosome here:
#microadriaticum: chr#_pilon
#kawagutii: HiC_scaffold_#
chromosome_num = "HiC_scaffold_1"

im = ax.matshow(
    clr.matrix(balance=True).fetch(chromosome_num),
    norm=norm,
    extent=(0,clr.chromsizes[chromosome_num], clr.chromsizes[chromosome_num], 0),
    
);
plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='raw counts');
ax.set_title('Chromosome ' + chromosome_num, y=1.08)
ax.set_ylabel('position, Mb')
format_ticks(ax)


region = (chromosome_num, 0, clr.chromsizes[chromosome_num])
extent=(0,clr.chromsizes[chromosome_num], clr.chromsizes[chromosome_num], 0),

cg = adaptive_coarsegrain(clr.matrix(balance=True).fetch(region),
                              clr.matrix(balance=False).fetch(region),
                              cutoff=3, max_levels=8)

cgi = interp_nan(cg)

f, axs = plt.subplots(
    figsize=(18,5),
    nrows=1,
    ncols=1,
    sharex=True, sharey=True)

im3 = axs.matshow(cgi, norm=norm, extent=extent)

#format_ticks(axs, rotate=False)
plt.gca().xaxis.tick_bottom()

plt.colorbar(im3, ax=axs, fraction=0.046, label='Normalized Contact Frequency')