### Request access to dataset hosted by Zenodo.
Create a zenodo account: https://zenodo.org/ you can login using your github account.\
Confirm your email attached to your zenodo account by clicking the link in the email sent to you by zenodo.\
Go to: https://zenodo.org/record/8277338. \
Request access to the data.\

### Install the following python packages:
```
pip install cooler
pip install cooltools
pip install numpy
pip install pandas
pip install matplotlib
```

### Very Important Google Drive Documents:

Miscellanious commands:
https://docs.google.com/document/d/1YIPwidYIA1nwAA2GbczWpTUEh0doGIe9/edit?usp=sharing&ouid=118151059404127538464&rtpof=true&sd=true

### Other Google Drive Documents:

Reading list:
https://docs.google.com/document/d/1HOtIcwAnNB_7rkMwDU-D4F7sM6BE9LPNg3ZfE19NNNA/edit?usp=sharing

Paper Outline:
https://docs.google.com/document/d/14jCeGCLy07exlhpo0kTJyaFfvMh-a48L/edit?usp=sharing&ouid=118151059404127538464&rtpof=true&sd=true

### Install MATLAB:
Install statistics and machine learning toolbox

### Other github repos used in this project:

https://github.com/lucasphilipp1/GEM
Description: Predicts 3D conformation from HiC data. Method is based on manifold learning, which contrasts CSynth. The 3D conformation take days to compute using this method.
Zhu, G., Deng, W., Hu, H., Ma, R., Zhang, S., Yang, J., ... & Zeng, J. (2018). Reconstructing spatial organizations of chromosomes through manifold learning. Nucleic acids research, 46(8), e50-e50.

https://github.com/lucasphilipp1/HIPPS-DIMES
Description: Estimates cell-cell conformational heterogeniety from HiC data.
Shi, G., & Thirumalai, D. (2023). A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants. Nature Communications, 14(1), 1150.

https://github.com/lucasphilipp1/3D_OrientationJ
Description: Image gradient based algorithm to calculate 3D orientatinon of DNA fibres from electron microscopy .tiff image z-stack. Unlike dragonfly, does not rely on segmentation to calculate ortientation.
Duclos, G., Adkins, R., Banerjee, D., Peterson, M. S., Varghese, M., Kolvin, I., ... & Dogic, Z. (2020). Topological structure and dynamics of three-dimensional active nematics. Science, 367(6482), 1120-1124.

https://github.com/diazale/dimension_reduction_workshop
Description: Workshop on t-SNE dimensionality reduction. Examples of how to implement t-SNE in python. Introductory level details to the theory are provided. Comparison with other dimensinoality reduction techniques. See in particular: toy_examples.ipynb.

https://github.com/ncbi/sra-tools
Description: Used for downloading transcriptomes in fastq format from the SRA (sequence read archive).

https://github.com/alexdobin/STAR/tree/master
Description: Use for aligning RNA-seq datasets to HiC assembly. See miscellanious commands document.

### Descriptions of code:

# automatic_CSynth.js
Description: Use to compute 3D conformations of multiple chromosomes, in succession, using CSynth. 
See miscellanious commands document.
To use, drag and drop onto CSynth tab in browser: https://csynth.github.io/csynth/csynth.html 

# chrom_length_order.py
Description: Finds the largest scaffolds or "chromosomes" from dinoflagellate HiC data. Suprisingly, scaffolds are not ordered perfectly by size.

# CSynth_FISH_Compare_Su_Cell_2020.m
Description: Use to assess accuracy of CSynth conformations and to optomize CSynth parameters, using a human cell line where HiC data and 3D FISH data (652 probes) exists.

# CSynth_FISH_Compare_Wang_Science_2016.m
Description: Older version of CSynth_FISH_Compare_Su_Cell_2020.m using a different FISH data set with fewer probes (30).

# dinoflagellate_Ps.m
Description: Possibly useless. Used to compute contact probability "P(s)" curves from dinoflagellate HiC data. Use the python cooltools package instead.

# fractal_equilbrium_load.m
Description: Used to compute contact probability "P(s)" curves from fractal and equilibrium globule conformations as a positive control for the calculation.

# GC_content.m
Description: Possibly useless. GC content uncertainty based on Ns or ambiguous bases is not calculated. Ns are prevalent in dinoflagellate HiC assemblies. Consider using: https://rdrr.io/cran/seqinr/man/GC.html instead. A better indicator of gene density would be to align RNA-seq reads to the assembly directly.  

# Cholesteric_HiC.m
Description: Generate cholesteric polymer model with extra chromosomal loops. Vary loop length, number of discs, cholesteric pitch, and other model parameters. Compute cholesteric HiC matrix.

# Cholesteric_HiC_test.m
Description: A copy of Cholesteric_HiC.m used to test implementation of intra-disc extra chromosomal loops, and variable loop lengths accomodated using cell arrays.

# constrained_RW_1D.m
Description: Function used to generate extrachromosomal loops. Random walks are generated separately for each of the spatial three dimensions.

# constrained_self_avoiding_RW_3D.m
Description: Function used to generate extrachromosomal loops. Combines random walks in each dimension to compute a 3D random walk. 

# contact_probability_xyz.m
Description: Function used to compute contact probability "P(s)" curves from a 3D structure. Calculation is different from computing contact probability curve from HiC matrix.

# Orientation_Order_Parameter.py
Description: Used to compute the correlation between two tangent vectors to the 3D chromosome structure separated by a certain primary sequence length, averaged over the entire chromosome.

# plot_P(s)_tree_of_life.py
Description: Used to compute the contact probability "P(s)" curves for many eukaryotic organisms across the tree of life.

# plot_P(s)_dinoflagellate.py
Description: Used to compute the contact probability "P(s)" curves for dinoflagellate HiC data testing affect of chromosome ordering (chrom_length_order.py) and differences in aseembly (Smic1.0 vs Smic1.1N).

# space_delimited_to_tab_delimited.m
Description: Possibly useless. Fixed an issue where .txt files written from cooler dump command were not properly delimited using tabs, as required by CSynth. An error in the bash command was fixed.

# straw.R
Description: Possibly useless. Was used to extract chromosome contacts from .hic format and output to .txt files. Made useless by hic2cool, see miscellanious commmands document. 

# GEM_plot.m
Description: Use to visualize 3D output of https://github.com/lucasphilipp1/GEM.

### Data Sources:
# IMR90 cells HiC data:
Rao, S. et al. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell 159, 1665–1680 (2014).

# IMR90 cells old FISH dataset:
Wang, S. et al. Spatial organization of chromatin domains and compartments in single chromosomes. Science 353, 598–602 (2016).

# IMR90 cells new FISH dataset:
Su, J. H., et al. (2020). Genome-scale imaging of the 3D organization and transcriptional activity of chromatin. Cell, 182(6), 1641-1659.

# GEM:
Zhu, Guangxiang, Wenxuan Deng, Hailin Hu, Rui Ma, Sai Zhang, Jinglin Yang, Jian Peng, Tommy Kaplan, and Jianyang Zeng. "Reconstructing spatial organizations of chromosomes through manifold learning." Nucleic acids research 46, no. 8 (2018): e50-e50.

# Fractal Globule & Equilibrium Globule Structures:
Lieberman-Aiden E., et al. Comprehensive mapping of long-range interactions reveals folding principles of the human genome. Science 2009 Oct 9;326(5950):289-93. PMID: 19815776. See: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199.
GEO Accession #: GSE18199

# HiC data for many eukaryotes:
From: Hoencamp, C., Dudchenko, O., Elbatsh, A. M., Brahmachari, S., Raaijmakers, J. A., van Schaik, T., ... & Rowland, B. D. (2021). 3D genomics across the tree of life reveals condensin II as a determinant of architecture type. Science, 372(6545), 984-989. \

GSM5182714 Tammar wallaby blood (Sample1096) \
GSM5182715 Domestic chicken splenic-derived lymphocytes (Sample4077) \
GSM5182716 Burmese python blood (Sample0838) \
GSM5182717 African clawed frog fibroblast (XTC) (Sample4078) \
GSM5182718 Red-bellied piranha muscle (Sample2033) \
GSM5182719 Brownbanded bamboo shark blood (Sample1575) \
GSM5182720 Arctic lamprey muscle (Sample3820) \
GSM5182721 Sea squirt tissue (Sample2441, Sample2444) \
GSM5182722 Sea squirt tissue (Sample2441) \
GSM5182723 European lancelet whole animal(s) (Sample1824) \
GSM5182724 Purple sea urchin muscle (Sample2581) \
GSM5182725 Yellow fever mosquito whole animal(s) (SAMN08028733, SAMN08028734, SAMN08028735) \
GSM5182726 Yellow fever mosquito whole animal(s) (SAMN08028734) \
GSM5182727 Southern house mosquito whole animal(s) (SAMN06546149) \
GSM5182728 Fruit fly whole animal(s) (Sample4079) \
GSM5182729 Tardigrade whole animal(s) (Sample3840) \
GSM5182730 Roundworm whole animal(s) (Sample4082) \
GSM5182731 Chinese liver fluke whole animal(s) (Sample2741) \
GSM5182732 California sea hare mix of CNS cells (Sample1819) \
GSM5182733 Moss animal whole animal(s) (Sample1845) \
GSM5182734 Acropora millepora stony coral branch fragments (Sample1884, Sample1992) \
GSM5182735 Sea gooseberry tissue (Sample1830) \
GSM5182736 Common mushroom fruiting body (Sample3828) \
GSM5182737 Baker's yeast cells (Sample4076) \
GSM5182738 Peanut leaves (Sample2455) \
GSM5182739 Bread wheat leaves (Sample4080, Sample4081) \
GSM5182740 Chinese muntjac blood (Sample1296) \
GSM5182741 Chinese muntjac fibroblasts (Sample2606) \
GSM5182742 Indian muntjac (Sample2604, Sample3923) \
GSM5182743 Indian muntjac fibroblasts (Sample3923) \
GSM5182744 Fruit fly PnM cell line (SAMN10242997) \

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102740
From Gibcus, J. H., Samejima, K., Goloborodko, A., Samejima, I., Naumova, N., Nuebler, J., ... & Dekker, J. (2018). A pathway for mitotic chromosome formation. Science, 359(6376), eaao6135. \
human (HeLa) cells, not synchronized: \
GSM2745898_17MAY11_HSC02RPAB-PE50_CCHiC-HeLa-NS-R1_23-05-2011_kittlere.5.10000.cool \
human (HeLa) cells, mitosis: \
GSM2745897_18FEB15_PE50_C66B1AC-A_Sample_HiC1-aka-CCHiC-HeLaS3CCL2p2-M-98.1000.cool

From: Marinov, G. K., Trevino, A. E., Xiang, T., Kundaje, A., Grossman, A. R., & Greenleaf, W. J. (2021). Transcription-dependent domain-scale three-dimensional genome organization in the dinoflagellate Breviolum minutum. Nature Genetics, 53(5), 613-617. \
GSM4658994	Dinoflagellate Breviolum minutum

From: Nand, A., Zhan, Y., Salazar, O. R., Aranda, M., Voolstra, C. R., & Dekker, J. (2021). Genetic and spatial organization of the unusual chromosomes of the dinoflagellate Symbiodinium microadriaticum. Nature Genetics, 53(5), 618-629. \
GSE152150	Dinoflagellate Symbiodinium microadriaticum coccoid Smic 1.1N

# 
https://www.science.org/doi/full/10.1126/science.aao6135
GEO Accession #: GSE102740

# 
GSM2801019	mC_S_minutum_25C

### TO DO LIST:

Fix total sequence length
	Vary: % of sequence in extrachromosomal loops
	 	 % in Inter-disc loops
		 % in Intra-disc loops
	Discs connected top to bottom consecutively
	Discs connected in random order
	Re-visiting discs
Add “thermal” noise

Orientation order parameter

Chiral Order Parameter
<img width="861" alt="image" src="https://github.com/lucasphilipp1/Dinoflagellate/assets/94249076/b35ade71-1599-4082-83b0-e19b5b968f32">

## High priority:
fullfile all Data imports in all MATLAB scripts: https://www.mathworks.com/help/matlab/ref/fullfile.html 

Cholestric_HiC code arrays need to be replaced with cell arrays (allow for variable size)
function to incorporate loops properly into primary sequence

Make a table of electron microscopy data vs dinoflagellate order (Lucas):
Discs perpendicular to long axis of chromosome?
Yes, No

HiC Data from Fugacium kawagutii? (Lucas)

## Medium priority:

diagonal relative error Su 2020

## Low priority:

Specifically, you will use computer simulations to:
- specifically test aspects of the cholesteric liquid crystal model for consistency with
dinoflagellate chromatin conformation capture (HiC) data
- assess uniqueness of dinoflagellate chromosome structure compared to other eukaryotes
- assess cell-to-cell chromosome conformational heterogeneity in dinoflagellates
- rigorously test the accuracy of inferred 3D structures
- map gene density, GC content, and other genomic annotations to 3D structures,
investigating structure-function relationships
  Figure 1. A sketch of how DNA is folded in the cholesteric liquid crystal chromosome model.

automatic cell counter software
 
- assess similarities/differences in chromosome structure across dinoflagellate species

RNA-seq data exists for both B minutum and S microadriaticum, much like https://www.nature.com/articles/s41467-023-35909-2 did with ATAC-seq reads we can map RNA-seq data to the CSynth structure by mapping the RNA-seq read count to the 3D location of the corresponding gene’s on the structure. we can PCA the CSynth monomers and make a plot of aggregate RNA-seq as a function of PC1 and PC2. We might see oscillations in the RNA-seq vs PC1 plot corresponding to inter-disc spacings for e.g.
no need for EU labelling


program to map RNA-seq, or any other sequences for that matter (like retrotransposon sequences from https://www.nature.com/articles/s41467-018-03724-9): onto the HiC assembly:
https://github.com/alexdobin/STAR
hypothesis being that retrotransposon and tandem repeats occupy the bulk of the chromosome because the “useful” genes are extruded out as extra-chromosomal loops

histogram of P(s) exponents

simulate missing sequence from cholesteric HiC matrix

retrotransposons

distribution of loop sizes.
noise to monomer positions

genomic distribution of 5hmC exists in B minutum: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2801019

so. i have been disturbed by the problem of cell to cell heterogeneity.
from emailing experimentalists in the HiC community, they seem not to trust a single “averaged” 3D structure inferred from HiC data when FISH data shows so much conformational variety from cell to cell. (see pic below)
i found this method: Distance Matrix to Ensemble of Structures or DIMES (https://www.nature.com/articles/s41467-023-36412-4) which uses information theory to calculate the minimally biased/maximum entropy joint distribution of positions of loci, P({xi}), whose average distance is equal to that given by the HiC data.
the approach is similar to this paper: https://doi.org/10.1073/pnas.1506257112 or https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.248101, but wolynes and zhang use information theory to calculate the minimally biased/maximum entropy potential energy function for molecular dynamics simulations. these involve long MD simulations in order for the polymer to sample the entire conformational energy landscape.
DIMES instead computes conformations directly from P({xi}), hence is quick, avoiding long MD simulations. it’s just a Lagrange multiplier calculation, followed by sampling from a distribution
DIMES cannot compute a “consensus” or “average” structure like CSynth, so all the other modelling we’ve done so far is still valid
DIMES is on github with good documentation
DIMES is written in python
i understand the information theory behind DIMES, i can explain it to the lab during a journal club
structural heterogeneity can be quantified using the average spread in the separation distribution P({xi}), averaged across all pairs of monomers, or the Q value (a measure of commonality of contacts amongst an ensemble of structures as discussed here https://doi.org/10.1073/pnas.1506257112).
I propose using:
human interphase HiC data as a - control, expecting high structural heterogeneity
human mitotic HiC data as a + control, expecting low structural heterogeneity
dinoflagellate chromosomes will be the experimental group, hopefully showing low structural heterogeneity given that we hypothesize permanently condensed chromosomes
if we can quantify the structural heterogeneity in dinoflagellates, and if it is low, we can:
support using CSynth as a “consensus” method, given the demonstrated low cell to cell variation (maybe there is only one structure to predict!)
provide evidence for permanent condensation in dinoflagellate chromosomes
Wolynes and Zhang found low structural heterogeneity at small genomic scales, but high structural heterogeneity at large genomic scales in humans
we could look for heterogeneity as a function of genomic separation using DIMES in dinoflagellates, expecting long range consensus to persist in the case of liquid crystal chromosomes
the DIMES paper also has a way of clustering conformations using t-sne. this method doesn’t require the length of the polymer to be the same, so we can compare structures from different organisms computed using CSynth. hopefully we see dinoflagellate chromosomes cluster together, but away from other eukaryotic chromosomes

keep in mind there are two optomization techniques
coarse grain structures if needed

t-sne DIMES paper Dnm equation
















