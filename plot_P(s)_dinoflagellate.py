"""
@author: lucasphilipp
a minor modification of:
https://cooltools.readthedocs.io/en/latest/notebooks/contacts_vs_distance.html
"""
import matplotlib.pyplot as plt
import pandas as pd
import cooler
import cooltools
import numpy as np

##############################################################################

import hicstraw

#Smic1.0
#microadriaticum: 94 chromosomes

# chrom_lengths_microadriaticum_Smic1_0=np.zeros(94+95)
# chrom_names_microadriaticum_Smic1_0=[]

# clr_microadriaticum_Smic1_0 = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files/Dinoflagellate HiC Data/GSE152150_Dino-HiC-D1plus-R1.smic1.0.mapq_30.1000.mcool::resolutions/5000')

# i=0
# for chrom in clr_microadriaticum_Smic1_0.chromnames:
#     chrom_lengths_microadriaticum_Smic1_0[i] = (clr_microadriaticum_Smic1_0.extent(chrom)[1]-clr_microadriaticum_Smic1_0.extent(chrom)[0])*1000
#     chrom_names_microadriaticum_Smic1_0.append(clr_microadriaticum_Smic1_0.chromnames[i])
#     i=i+1

# chrom_length_ind_microadriaticum_Smic1_0 = np.argsort(chrom_lengths_microadriaticum_Smic1_0)
# chrom_length_ind_microadriaticum_Smic1_0 = chrom_length_ind_microadriaticum_Smic1_0[::-1]

# chrom_names_microadriaticum_Smic1_0 = [chrom_names_microadriaticum_Smic1_0[i] for i in chrom_length_ind_microadriaticum_Smic1_0]

# chrom_lengths_microadriaticum_Smic1_0 = [chrom_lengths_microadriaticum_Smic1_0[i] for i in chrom_length_ind_microadriaticum_Smic1_0]

#SMIC1.1N
#microadriaticum: 94 chromosomes

chrom_lengths_microadriaticum_Smic1_1N=np.zeros(100)
chrom_lengths_microadriaticum_Smic1_1N=np.zeros(94)
chrom_names_microadriaticum_Smic1_1N=[]

clr_microadriaticum_Smic1_1N = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/GSE152150_HiC-Dplus.smic1.1N.mapq_30.1000.mcool::/resolutions/1000')

i=0
for chrom in clr_microadriaticum_Smic1_1N.chromnames:
    chrom_lengths_microadriaticum_Smic1_1N[i] = (clr_microadriaticum_Smic1_1N.extent(chrom)[1]-clr_microadriaticum_Smic1_1N.extent(chrom)[0])*1000
    chrom_names_microadriaticum_Smic1_1N.append(clr_microadriaticum_Smic1_1N.chromnames[i])
    i=i+1

chrom_length_ind_microadriaticum_Smic1_1N = np.argsort(chrom_lengths_microadriaticum_Smic1_1N)
chrom_length_ind_microadriaticum_Smic1_1N = chrom_length_ind_microadriaticum_Smic1_1N[::-1]

chrom_names_microadriaticum_Smic1_1N = [chrom_names_microadriaticum_Smic1_1N[i] for i in chrom_length_ind_microadriaticum_Smic1_1N]

chrom_lengths_microadriaticum_Smic1_1N = [chrom_lengths_microadriaticum_Smic1_1N[i] for i in chrom_length_ind_microadriaticum_Smic1_1N]

#minutum: 91 pseudochromosomes > 500kb
#longest being ∼11megabases (Mb) in size, with a median length of 6.7Mb

chrom_lengths_minutum = np.zeros(16003)
chrom_names_minutum = []

hic = hicstraw.HiCFile('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/GSM4658994_L142+L533+L534.inter.hic')

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

print(len([i for i in chrom_lengths_minutum if i > 500000]))
print([chrom_names_minutum[chrom_lengths_minutum.index(i)] for i in chrom_lengths_minutum if i > 500000])
print(chrom_names_minutum[91])

##############################################################################

resolution_microadriaticum  = 5000
resolution_minutum  = 5000

clr_microadriaticum = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/GSE152150_HiC-Dplus.smic1.1N.mapq_30.1000.mcool::/resolutions/5000')
#clr_microadriaticum = cooler.Cooler('/Users/lucasphilipp/Downloads/GSE152150_Dino-HiC-cD3-DSG-TPL-DRB.smic1.1N.mapq_30.1000.mcool::/resolutions/5000')

#clr_microadriaticum = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/Symbiodinium_kawagutii_5000.cool')

clr_minutum = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/GSM4658994_L142+L533+L534.inter_first_100_scaffolds.mcool::/resolutions/5000')

df_microadriaticum = pd.DataFrame({'chrom': clr_microadriaticum.chromnames,
'start': 0,
'end': clr_microadriaticum.chromsizes.values,
'name': clr_microadriaticum.chromnames}
)

df_minutum = pd.DataFrame({'chrom': clr_minutum.chromnames,
'start': 0,
'end': clr_minutum.chromsizes.values,
'name': clr_minutum.chromnames}
)

#f_microadriaticum=df_microadriaticum.loc[chrom_length_ind_microadriaticum_Smic1_1N[0:94]]
#df_microadriaticum=df_microadriaticum.loc[chrom_length_ind_microadriaticum_Smic1_1N[0:99]]
#df_minutum = df_minutum.loc[chrom_length_ind_minutum[0:91]]

# cvd == contacts-vs-distance
cvd_smooth_microadriaticum = cooltools.expected_cis(
clr=clr_microadriaticum,
view_df=df_microadriaticum,
smooth=True,
aggregate_smoothed=True,
nproc=8 #if you do not have multiple cores available, set to 1
)

cvd_smooth_minutum = cooltools.expected_cis(
clr=clr_minutum,
view_df=df_minutum,
smooth=True,
aggregate_smoothed=True,
nproc=8 #if you do not have multiple cores available, set to 1
)

#view output of P(s) calculation
#print(cvd_smooth_microadriaticum.columns)
#print(cvd_smooth_minutum.columns)

#columns are chr, chr, s, Ps_smoothed
#print(cvd_smooth_microadriaticum.values[:, [0,1,2,8]])
#print(cvd_smooth_minutum.values[:, [0,1,2,8]])

#P(s) curve for raw counts: count.avg, P(s) curve for normalized counts balanced.avg
#pd.DataFrame(cvd_smooth_microadriaticum.values[:, [0,1,2,8]]).to_csv('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Ps_symbiodinium_microadriaticum_coccoid/GSE152150_Dino-HiC-cD1plus-R1.smic1.1N.mapq_30.1000.mcool_P(s).txt', sep='\t')
#pd.DataFrame(cvd_smooth_minutum.values[:, [0,1,2,8]]).to_csv('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Ps_bminutum/GSM4658994_L142+L533+L534.inter.mcool_P(s).txt', sep='\t')

cvd_smooth_microadriaticum['s_bp'] = cvd_smooth_microadriaticum['dist'] * resolution_microadriaticum
cvd_smooth_microadriaticum['balanced.avg.smoothed'].loc[cvd_smooth_microadriaticum['dist'] < 2] = np.nan

cvd_smooth_minutum['s_bp'] = cvd_smooth_minutum['dist'] * resolution_minutum
cvd_smooth_minutum['balanced.avg.smoothed'].loc[cvd_smooth_minutum['dist'] < 2] = np.nan

# for region in df_microadriaticum['name']:
#     f, ax = plt.subplots(1,1)
#     ax.loglog(
#     cvd_smooth_microadriaticum['s_bp'].loc[cvd_smooth_microadriaticum['region1']==region],
#     cvd_smooth_microadriaticum['balanced.avg'].loc[cvd_smooth_microadriaticum['region1']==region],
#     cvd_smooth_microadriaticum['s_bp'].loc[cvd_smooth_microadriaticum['region1']==region],
#     cvd_smooth_microadriaticum['balanced.avg.smoothed'].loc[cvd_smooth_microadriaticum['region1']==region]
#     )
#     ax.set(
#     xlabel='s, bp',
#     ylabel='P(s)')
#     ax.set_aspect(1.0)
#     ax.grid(lw=0.5)
#     #plt.savefig("/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Ps_symbiodinium_microadriaticum_coccoid/Plots_{0}".format(region))
#     plt.show()
    
# for region in df_minutum['name']:
#     f, ax = plt.subplots(1,1)
#     ax.loglog(
#     cvd_smooth_minutum['s_bp'].loc[cvd_smooth_minutum['region1']==region],
#     cvd_smooth_minutum['balanced.avg'].loc[cvd_smooth_minutum['region1']==region],
#     cvd_smooth_minutum['s_bp'].loc[cvd_smooth_minutum['region1']==region],
#     cvd_smooth_minutum['balanced.avg.smoothed'].loc[cvd_smooth_minutum['region1']==region]
#     )
#     ax.set(
#     xlabel='s, bp',
#     ylabel='P(s)')
#     ax.set_aspect(1.0)
#     ax.grid(lw=0.5)
#     #plt.savefig("/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Ps_bminutum/Plots_{0}".format(region))
#     plt.show()

x_line = np.linspace(10**4, 10**7, 100)
y_line = (x_line**-0.2)/110
y_line2 = (x_line**-0.4)/9
    
f, axs = plt.subplots(2, figsize=(15, 15))
ax=axs[0]
for region in df_microadriaticum['name']:
    ax.loglog(
    cvd_smooth_microadriaticum['s_bp'].loc[cvd_smooth_microadriaticum['region1']==region],
    cvd_smooth_microadriaticum['balanced.avg.smoothed'].loc[cvd_smooth_microadriaticum['region1']==region], alpha=0.3, color='blue', label = 'microadriaticum'
    )
    #ax.loglog([x_line[0], x_line[-1]], [y_line[0], y_line[-1]], 'r--')
    ax.loglog([x_line[0], x_line[-1]], [y_line2[0], y_line2[-1]], 'k--')

    
# for region in df_minutum['name']:
#     ax.loglog(
#     cvd_smooth_minutum['s_bp'].loc[cvd_smooth_minutum['region1']==region],
#     cvd_smooth_minutum['balanced.avg.smoothed'].loc[cvd_smooth_minutum['region1']==region], alpha=0.3, color='orange', label = 'minutum'
#     )
#     ax.loglog([x_line[0], x_line[-1]], [y_line2[0], y_line2[-1]], 'k--')

# handles, labels = ax.get_legend_handles_labels()
# display = (0,94)
# ax.legend([handle for i,handle in enumerate(handles) if i in display],
#       [label for i,label in enumerate(labels) if i in display], loc = 'best')
ax.set(
xlabel='s, bp',
ylabel='P(s)')
ax.set_ylim(10**-6,10**0)
ax.set_aspect(0.55)
ax.grid(lw=0.5)    
plt.xlim( (10**4,2*10**7) )
#plt.ylim( (10**-8,10**0) )

grad_all = np.zeros(1);

ax=axs[1]
# for region in df_microadriaticum['name']:
#     y=cvd_smooth_microadriaticum['balanced.avg.smoothed'].loc[cvd_smooth_microadriaticum['region1']==region]
#     x=cvd_smooth_microadriaticum['s_bp'].loc[cvd_smooth_microadriaticum['region1']==region]
#     logy=np.log(y)
#     logx=np.log(x)
#     grad=np.gradient(logy,logx)
    
#     grad_all = np.append(grad_all, grad[0:max(max(np.where(x.to_numpy()<1*10**6)))])

#     # Calculate derivative in log-log space
#     ax.semilogx(
#     cvd_smooth_microadriaticum['s_bp'].loc[cvd_smooth_microadriaticum['region1']==region],
#     grad, alpha=0.3, color='blue', label = 'microadriaticum')

   
for region in df_minutum['name']:
    y=cvd_smooth_minutum['balanced.avg.smoothed'].loc[cvd_smooth_minutum['region1']==region]
    x=cvd_smooth_minutum['s_bp'].loc[cvd_smooth_minutum['region1']==region]
    logy=np.log(y)
    logx=np.log(x)
    grad=np.gradient(logy,logx)
    # Calculate derivative in log-log space
    ax.semilogx(
    cvd_smooth_minutum['s_bp'].loc[cvd_smooth_minutum['region1']==region],
    grad, alpha=0.3, color='orange', label = 'minutum')
    
# handles, labels = ax.get_legend_handles_labels()
# display = (0,94)
# ax.legend([handle for i,handle in enumerate(handles) if i in display],
#       [label for i,label in enumerate(labels) if i in display], loc = 'best')
ax.set(
xlabel='s, bp',
ylabel='slope')
ax.set_aspect(1.04)
ax.grid(lw=0.5)

plt.rcParams.update({'font.size': 25})
#plt.xlim( (7*10**3,2.75*10**7) )
plt.xlim( (10**4,2*10**7) )
plt.ylim( (-2,1) )
plt.show()

# grad_all = grad_all[~np.isnan(grad_all)]
# n, bins, patches = plt.hist(grad_all)
# plt.show()
# mode_index = n.argmax()    
# print('the mode:'+ str((bins[mode_index] + bins[mode_index+1])/2))


    
    