"""
@author: lucasphilipp
a modification of:
https://cooltools.readthedocs.io/en/latest/notebooks/contacts_vs_distance.html
"""

#columns are chr, chr, s, Ps_smoothed
#print(cvd_smooth.values[:, [0,1,2,8]])
#P(s) curve for raw counts: count.avg, P(s) curve for normalized counts balanced.avg

#choose a file path to save the P(s) curve as a .csv file
#pd.DataFrame(cvd_smooth.values[:, [0,1,2,8]]).to_csv(filepath + '_' + organism + '_P(s).csv', sep='\t')

import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import cooler
import cooltools
import numpy as np

class organism_HiC:
    def __init__(self, filepath: str, resolution: int=5000)->None:
        self.filepath = filepath
        self.resolution = resolution
        self.clr = cooler.Cooler(filepath)
        self.df = pd.DataFrame({'chrom': self.clr.chromnames,
        'start': 0,
        'end': self.clr.chromsizes.values,
        'name': self.clr.chromnames}
        )

        self.c = cooltools.expected_cis(
        clr=self.clr,
        view_df=self.df,
        smooth=True,
        aggregate_smoothed=True,
        nproc=8 #if you do not have multiple cores available, set to 1
        )
        
        self.c['s_bp'] = self.c['dist'] * self.resolution
        self.c['balanced.avg.smoothed'].loc[self.c['dist'] < 2] = np.nan
        self.c['balanced.avg.smoothed.agg'].loc[self.c['dist'] < 2] = np.nan
    
def utility_display(org_list):
    disp_arr = np.zeros(len(org_list))
    for i in range(len(org_list)-1):
        disp_arr[i+1] = len(org_list[i].clr.chromnames)
    #return the cum_sum
    return np.cumsum(disp_arr)


def batch_plot(figure_obj,ax_tuple,org_list,label_list):
    ax = ax_tuple
    for i in range(len(org_list)):
        for region in org_list[i].df['name']:
            ax[0].loglog(
            org_list[i].c['s_bp'].loc[org_list[i].c['region1']==region],
            org_list[i].c['balanced.avg.smoothed.agg'].loc[org_list[i].c['region1']==region],
            #alpha=0.3, color = cm.rainbow(i/len(org_list)), label = label_list[i]
            alpha=0.3, color = cm.rainbow((30-i)/30), label = label_list[i]
            )
            y=org_list[i].c['balanced.avg.smoothed.agg'].loc[org_list[i].c['region1']==region]
            x=org_list[i].c['s_bp'].loc[org_list[i].c['region1']==region]
            logy=np.log(y)
            logx=np.log(x)
            grad=np.gradient(logy,logx)
            # Calculate derivative in log-log space
            ax[1].semilogx(
            org_list[i].c['s_bp'].loc[org_list[i].c['region1']==region],
            #grad, alpha=0.3, color = cm.rainbow(i/len(org_list)), label = label_list[i]
            grad, alpha=0.3, color = cm.rainbow((30-i)/30), label = label_list[i]
            )
    ax[0].set(
    xlabel='s, bp',
    ylabel='P(s)')
    ax[0].grid(lw=0.5)
    ax[0].set_aspect(1.0)

    plt.xlim( (10**3,10**8) )
    plt.ylim( (10**-6,10**-1) )
    
    handles, labels = ax[0].get_legend_handles_labels()
    
    display = utility_display(org_list)
    
    ax[1].legend([handle for i,handle in enumerate(handles) if i in display],
          [label for i,label in enumerate(labels) if i in display], loc='upper center', bbox_to_anchor=(1.5, 1),
                    fancybox=True, shadow=True, ncol=1)
    ax[1].set(
    xlabel='s, bp',
    ylabel='slope')
    ax[1].grid(lw=0.5)
    ax[1].set_aspect(1.0)

    plt.xlim( (10**3,10**8) )
    plt.ylim( (-3,2) )
        
    plt.show()

def chr_variation_per_organism_plot(org_list,label_list):
    #this loop displays the raw P(s) curves and smoothed P(s) curves for each chromosome
    for i in range(len(org_list)):
        for region in org_list[i].df['name']:
            f, ax = plt.subplots(1,1)
            ax.loglog(
            org_list[i].c['s_bp'].loc[org_list[i].c['region1']==region],
            org_list[i].c['balanced.avg'].loc[org_list[i].c['region1']==region],
            org_list[i].c['s_bp'].loc[org_list[i].c['region1']==region],
            org_list[i].c['balanced.avg.smoothed'].loc[org_list[i].c['region1']==region],
            label = label_list[i]
            )
            ax.set(
            xlabel='s, bp',
            ylabel='P(s)') 
            ax.set_aspect(1.0)
            ax.grid(lw=0.5)
            plt.show()

org_list = []

org1 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Acropora millepora stony coral branch fragments/GSM5182734_amil_sf_1.1_HiC.mcool::resolutions/5000')
org_list.append(org1)
print(1)

org2 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/African clawed frog fibroblast/GSM5182717_Xla.v91.mcool::resolutions/5000')
org_list.append(org2)
print(2)

org3 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Arctic lamprey muscle/GSM5182720_LetJap1.0_HiC.mcool::resolutions/5000')
org_list.append(org3)
print(3)

org4 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Bakers yeast/GSM5182737_sacCer3.mcool::resolutions/5000')
org_list.append(org4)
print(4)

org5 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Brownbanded bamboo shark blood/GSM5182719_Cpunctatum_v1.0_HiC.mcool::resolutions/5000')
org_list.append(org5)
print(5)

org6 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Burmese python blood/GSM5182716_Python_molurus_bivittatus-5.0.2_HiC.mcool::resolutions/5000')
org_list.append(org6)
print(6)

org7 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/California sea hare mix of CNS cells/GSM5182732_AplCal3.0_HiC.mcool::resolutions/5000')
org_list.append(org7)
print(7)

org8 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Chinese liver fluke whole animal(s)/GSM5182731_ASM360417v1_HiC.mcool::resolutions/5000')
org_list.append(org8)
print(8)

org9 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Chinese muntjac/GSM5182741_M.reevesi_fibroblast_vs_CIJ_HiC_reordered.mcool::resolutions/5000')
org_list.append(org9)
print(9)

org10 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Common mushroom fruiting body/GSM5182736_Agabi_varbisH97_2_HiC.mcool::resolutions/5000')
org_list.append(org10)
print(10)

org11 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Domestic chicken splenic-derived lymphocytes/GSM5182715_GRCg6a.mcool::resolutions/5000')
org_list.append(org11)
print(11)

org12 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/European lancelet whole animal(s)/GSM5182723_Bl71nemr_HiC.mcool::resolutions/5000')
org_list.append(org12)
print(12)

org13 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Fruit fly whole animal(s)/GSM5182728_Release_6_plus_ISO1_MT_merged.mcool::resolutions/5000')
org_list.append(org13)
print(13)

org14 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Indian muntjac/GSM5182743_M.muntjak_fibroblast_vs_CMJ_HiC_female.mcool::resolutions/5000')
org_list.append(org14)
print(14)

org15 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Moss animal whole animal(s)/GSM5182733_cmucedo.flye2_HiC.mcool::resolutions/5000')
org_list.append(org15)
print(15)

org16 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Peanut leaves/GSM5182738_arahy.Tifrunner.gnm2.J5K5.mcool::resolutions/5000')
org_list.append(org16)
print(16)

org17 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Purple sea urchin muscle/GSM5182724_Spur_5.0.contigs_purged_HiC.mcool::resolutions/5000')
org_list.append(org7)
print(17)

org18 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Red-bellied piranha muscle/GSM5182718_Pygocentrus_nattereri-1.0.2_HiC.mcool::resolutions/5000')
org_list.append(org18)
print(18)

org19 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Roundworm whole animal(s)/GSM5182730_WBcel235.mcool::resolutions/5000')
org_list.append(org19)
print(19)

org20 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Sea gooseberry tissue/GSM5182735_P.bachei_draft_genome_v.1.1_HiC.mcool::resolutions/5000')
org_list.append(org20)
print(20)

org21 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Sea squirt tissue/GSM5182721_Ciona_intestinalis_HT_Hoya_T-line_assembly_2019_HiC.mcool::resolutions/5000')
org_list.append(org21)
print(21)

org22 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Southern house mosquito whole animal(s)/GSM5182727_CpipJ3.mcool::resolutions/5000')
org_list.append(org22)
print(22)

org23 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Tammar wallaby blood/GSM5182714_me-1k.mcool::resolutions/5000')
org_list.append(org23)
print(23)

org24 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Tardigrade whole animal(s)/GSM5182729_nHd_3.1_HiC.mcool::resolutions/5000')
org_list.append(org24)
print(24)

org25 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Yellow fever mosquito whole animal(s)/GSM5182725_AaegL5.0.mcool::resolutions/5000')
org_list.append(org25)
print(25)

org27 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Human/GSM2745898_17MAY11_HSC02RPAB-PE50_CCHiC-HeLa-NS-R1_23-05-2011_kittlere.5.10000.cool')
org_list.append(org27)
print(27)

org28 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Human/GSM2745897_18FEB15_PE50_C66B1AC-A_Sample_HiC1-aka-CCHiC-HeLaS3CCL2p2-M-98.1000.cool')
org_list.append(org28)
print(28)

org29 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Dinoflagellate Breviolum minutum/GSM4658994_L142+L533+L534.inter_first_100_scaffolds.mcool::/resolutions/5000')
org_list.append(org29)
print(29)

org30 = organism_HiC('/Users/lucasphilipp/Desktop/Research/GitHub/HiC Across The Tree of Life/Dinoflagellate Symbiodinium Microadriaticum/GSE152150_Dino-HiC-cD1plus-R1.smic1.1N.mapq_30.1000.mcool::/resolutions/5000')
org_list.append(org30)
print(30)

label_list=[]

label1 = 'Acropora millepora stony coral branch fragments'
label2 = 'African clawed frog fibroblast'
label3 = 'Arctic lamprey muscle'
label4 = 'Bakers yeast'
label5 = 'Brownbanded bamboo shark blood'
label6 = 'Burmese python blood'
label7 = 'California sea hare mix of CNS cells'
label8 = 'Chinese liver fluke whole animal(s)'
label9 = 'Chinese muntjac'
label10 = 'Common mushroom fruiting body'
label11 = 'Domestic chicken splenic-derived lymphocytes'
label12 = 'European lancelet whole animal(s)'
label13 = 'Fruit fly whole animal(s)'
label14 = 'Indian muntjac'
label15 = 'Moss animal whole animal(s)'
label16 = 'Peanut leaves'
label17 = 'Purple sea urchin muscle'
label18 = 'Red-bellied piranha muscle'
label19 = 'Roundworm whole animal(s)'
label20 = 'Sea gooseberry tissue'
label21 = 'Sea squirt tissue'
label22 = 'Southern house mosquito whole animal(s)'
label23 = 'Tammar wallaby blood'
label24 = 'Tardigrade whole animal(s)'
label25 = 'Yellow fever mosquito whole animal(s)'
label27 = 'Human HeLa: not synchronized'
label28 = 'Human HeLa: mitosis'
label29 = 'Dinoflagellate: Breviolum minutum'
label30 = 'Dinoflagellate: Symbiodinium microadriaticum'

label_list.append(label1)
label_list.append(label2)
label_list.append(label3)
label_list.append(label4)
label_list.append(label5)
label_list.append(label6)
label_list.append(label7)
label_list.append(label8)
label_list.append(label9)
label_list.append(label10)
label_list.append(label11)
label_list.append(label12)
label_list.append(label13)
label_list.append(label14)
label_list.append(label15)
label_list.append(label16)
label_list.append(label17)
label_list.append(label18)
label_list.append(label19)
label_list.append(label20)
label_list.append(label21)
label_list.append(label22)
label_list.append(label23)
label_list.append(label24)
label_list.append(label25)
label_list.append(label27)
label_list.append(label28)
label_list.append(label29)
label_list.append(label30)

# color1 = 'black'
# color2 = 'gray'
# color3 = 'rosybrown'
# color4 = 'brown'
# color5 = 'red'
# color6 = 'darkorange'
# color7 = 'orange'
# color8 = 'gold'
# color9 = 'darkkhaki'
# color10 ='olive'
# color11 ='chartreuse'
# color12 ='palegreen'
# color13 ='darkgreen'
# color14 ='springgreen'
# color15 ='aquamarine'
# color16 ='teal'
# color17 ='lightblue'
# color18 ='deepskyblue'
# color19 ='royalblue'
# color20 ='mediumpurple'
# color21 ='indigo'
# color22 ='violet'
# color23 ='hotpink'
# color24 ='crimson'
# color25 ='lightpink'

#f, (ax, ax2) = plt.subplots(1,2, figsize=(12, 12))
#batch_plot(f, (ax, ax2), org_list, label_list)
chr_variation_per_organism_plot(org_list, label_list)
