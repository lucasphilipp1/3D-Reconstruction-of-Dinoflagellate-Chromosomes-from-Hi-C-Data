import hicstraw
import pandas as pd

#hic_file = 'Symbiodinium_kawagutii.hic'
hic_file = 'GSE113256_mega.hic'

data_type = 'observed' # (previous default / "main" data) or 'oe' (observed/expected)
normalization = 'KR'  # , VC, VC_SQRT, KR, SCALE, etc.

resolution = 25000

hic = hicstraw.HiCFile(hic_file)
print(hic.getGenomeID())

assert resolution in hic.getResolutions(), \
    f"{resolution} is not part of the possible resolutions {','.join(hic.getResolutions())}"

chrom_sizes = pd.Series({chrom.name: chrom.length for chrom in hic.getChromosomes() if chrom.name != "All"})
print(chrom_sizes.to_string())

#replace with number of HIC_SCAFFOLDS
for i in range(1):
        chrom = chrom_sizes.index[i]
        print(chrom)
        result = hicstraw.straw(data_type, normalization, hic_file, 'NC_035107.1', 'NC_035107.1', 'BP', resolution)
        #result = hicstraw.straw(data_type, normalization, hic_file, chrom, chrom, 'BP', resolution)
        
        #with open('s_kawagutii_V3_'+chrom+'.txt','w') as file:
        with open('mosquito'+'_NC_035107.1_25000.txt','w') as file:
            for k in range(len(result)):
                start1 = result[k].binX
                start2 = result[k].binY
                value = result[k].counts
                file.write(str(start1) + '\t' + str(start2) + '\t' + str(value) + '\n')

