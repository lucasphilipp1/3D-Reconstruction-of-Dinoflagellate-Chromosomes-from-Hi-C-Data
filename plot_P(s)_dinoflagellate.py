"""
a minor modification of:
https://cooltools.readthedocs.io/en/latest/notebooks/contacts_vs_distance.html

Fig 1C and Fig S7 in the paper
"""
import matplotlib.pyplot as plt
import pandas as pd
import cooler
import cooltools
import numpy as np

##############################################################################

resolution_microadriaticum  = 5000
resolution_kawagutii  = 5000
resolution_minutum  = 5000

#clr_microadriaticum = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/GSE152150_HiC-Dplus.smic1.1N.mapq_30.1000.mcool::/resolutions/5000')
#clr_kawagutii = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/Symbiodinium_kawagutii_5000.cool')
clr_minutum = cooler.Cooler('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/GSM4658994_L142+L533+L534.inter_first_100_scaffolds.mcool::/resolutions/5000')

# df_microadriaticum = pd.DataFrame({'chrom': clr_microadriaticum.chromnames,
# 'start': 0,
# 'end': clr_microadriaticum.chromsizes.values,
# 'name': clr_microadriaticum.chromnames}
# )

# df_kawagutii = pd.DataFrame({'chrom': clr_kawagutii.chromnames,
# 'start': 0,
# 'end': clr_kawagutii.chromsizes.values,
# 'name': clr_kawagutii.chromnames}
# )

df_minutum = pd.DataFrame({'chrom': clr_minutum.chromnames,
'start': 0,
'end': clr_minutum.chromsizes.values,
'name': clr_minutum.chromnames}
)

#df_microadriaticum=df_microadriaticum.loc[chrom_length_ind_microadriaticum_Smic1_1N[0:94]]
#df_kawagutii=df_kawagutii.loc[chrom_length_ind_microadriaticum_Smic1_1N[0:99]]
#df_minutum = df_minutum.loc[chrom_length_ind_minutum[0:91]]

# # cvd == contacts-vs-distance
# cvd_smooth_microadriaticum = cooltools.expected_cis(
# clr=clr_microadriaticum,
# view_df=df_microadriaticum,
# smooth=True,
# aggregate_smoothed=True,
# nproc=8 #if you do not have multiple cores available, set to 1
# )

# # cvd == contacts-vs-distance
# cvd_smooth_kawagutii = cooltools.expected_cis(
# clr=clr_kawagutii,
# view_df=df_kawagutii,
# smooth=True,
# aggregate_smoothed=True,
# nproc=8 #if you do not have multiple cores available, set to 1
# )

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

#ignore main diagonal and first off-diagonal due to self-ligation and unligated artefacts
# cvd_smooth_microadriaticum['s_bp'] = cvd_smooth_microadriaticum['dist'] * resolution_microadriaticum
# cvd_smooth_microadriaticum['balanced.avg.smoothed'].loc[cvd_smooth_microadriaticum['dist'] < 2] = np.nan

# cvd_smooth_microadriaticum['s_bp'] = cvd_smooth_microadriaticum['dist'] * resolution_microadriaticum
# cvd_smooth_microadriaticum['balanced.avg.smoothed'].loc[cvd_smooth_microadriaticum['dist'] < 2] = np.nan

cvd_smooth_minutum['s_bp'] = cvd_smooth_minutum['dist'] * resolution_minutum
cvd_smooth_minutum['balanced.avg.smoothed'].loc[cvd_smooth_minutum['dist'] < 2] = np.nan

#add -0.5 and -0.2 exponents to plot

x_line = np.linspace(10**4, 10**7, 100)
y_line = (x_line**-0.2)/140
y_line2 = (x_line**-0.5)/3.5
    
f, axs = plt.subplots(2, figsize=(15, 15))
ax=axs[0]
# for region in df_microadriaticum['name']:
#     ax.loglog(
#     cvd_smooth_microadriaticum['s_bp'].loc[cvd_smooth_microadriaticum['region1']==region],
#     cvd_smooth_microadriaticum['balanced.avg.smoothed'].loc[cvd_smooth_microadriaticum['region1']==region], alpha=0.3, color='blue', label = 'microadriaticum'
#     )
#     #ax.loglog([x_line[0], x_line[-1]], [y_line[0], y_line[-1]], 'r--')
#     ax.loglog([x_line[0], x_line[-1]], [y_line2[0], y_line2[-1]], 'k--')

    
for region in df_minutum['name']:
    ax.loglog(
    cvd_smooth_minutum['s_bp'].loc[cvd_smooth_minutum['region1']==region],
    cvd_smooth_minutum['balanced.avg.smoothed'].loc[cvd_smooth_minutum['region1']==region], alpha=0.3, color='orange', label = 'minutum'
    )
    ax.loglog([x_line[0], x_line[-1]], [y_line2[0], y_line2[-1]], 'k--')
    ax.loglog([x_line[0], x_line[-1]], [y_line[0], y_line[-1]], 'r--')

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

    
    