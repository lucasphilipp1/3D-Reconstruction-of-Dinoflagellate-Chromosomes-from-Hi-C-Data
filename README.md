ALL STATIC DATA SETS STORED HERE: 

### Request access to dataset hosted by Zenodo.
Create a zenodo account: https://zenodo.org/ you can login using your github account.\
Confirm your email attached to your zenodo account by clicking the link in the email sent to you by zenodo.\
Go to: https://zenodo.org/record/8277338. \
Request access to the data.\
Once approved to access the dataset, choose your favourite organism and download the appropriate .zip file.\
Open the folder. You should see a .hic file in the folder.\
If your folder doesn't have a .hic file don't worry, download another organism from zenodo.

### Install the following python packages:
```
pip install cooler
pip install cooltools
pip install numpy
pip install pandas
pip install matplotlib
```

### Install MATLAB

Other github repos used in this project:

https://github.com/lucasphilipp1/GEM
Description:

https://github.com/lucasphilipp1/HIPPS-DIMES
Description:

https://github.com/lucasphilipp1/3D_OrientationJ
Description:

### Descriptions of code:

automatic_CSynth.js
Description:

chrom_length_order.py
Description:

CSynth_FISH_Compare_Su_Cell_2020.m
Description:

CSynth_FISH_Compare_Wang_Science_2016.m
Description:

dinoflagellate_Ps.m
Description:

fractal_equilbrium_load.m
Description:

GC_content.m
Description:

Cholesteric_HiC.m
Description:

constrained_RW_1D.m
Description:

constrained_self_avoiding_RW_3D.m
Description:

contact_probability_xyz.m
Description:

Cholesteric_HiC_test.m
Description:

Orientation_Order_Parameter.py
Description:

plot_P(s)_workshop.py
Description:

plot_P(s).py
Description:

space_delimited_to_tab_delimited.m
Description:

straw.R
Description:

GEM_plot.m
Description:

### Data Sources:
IMR90 cells HiC data:
Rao, S. et al. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell 159, 1665–1680 (2014).

IMR90 cells FISH data:
Wang, S. et al. Spatial organization of chromatin domains and compartments in single chromosomes. Science 353, 598–602 (2016).

GEM:
Zhu, Guangxiang, Wenxuan Deng, Hailin Hu, Rui Ma, Sai Zhang, Jinglin Yang, Jian Peng, Tommy Kaplan, and Jianyang Zeng. "Reconstructing spatial organizations of chromosomes through manifold learning." Nucleic acids research 46, no. 8 (2018): e50-e50.

Fractal Globule & Equilibrium Globule Structures:
Lieberman-Aiden E*, van Berkum NL*, Williams L, Imakaev M et al. Comprehensive mapping of long-range interactions reveals folding principles of the human genome. Science 2009 Oct 9;326(5950):289-93. PMID: 19815776. See: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199.



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
Cholestric_HiC code arrays need to be replaced with cell arrays (allow for variable size)
function to incorporate loops properly into primary sequence

## Medium priority:

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
 
- assess similarities/differences in chromosome structure across dinoflagellate species

RNA-seq data exists for both B minutum and S microadriaticum, much like https://www.nature.com/articles/s41467-023-35909-2 did with ATAC-seq reads we can map RNA-seq data to the CSynth structure by mapping the RNA-seq read count to the 3D location of the corresponding gene’s on the structure. we can PCA the CSynth monomers and make a plot of aggregate RNA-seq as a function of PC1 and PC2. We might see oscillations in the RNA-seq vs PC1 plot corresponding to inter-disc spacings for e.g.
no need for EU labelling (edited) 
:grinning:
1

7:59
program to map RNA-seq, or any other sequences for that matter (like retrotransposon sequences from https://www.nature.com/articles/s41467-018-03724-9): onto the HiC assembly:
https://github.com/alexdobin/STAR
hypothesis being that retrotransposon and tandem repeats occupy the bulk of the chromosome because the “useful” genes are extruded out as extra-chromosomal loops

histogram of P(s) exponents

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











