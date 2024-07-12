# Zenodo

# Viewing HiC Data:
.cool & .mcool format:
https://resgen.io/
account creation is required to download a dataset
.hic format:
https://www.aidenlab.org/juicebox/

# Very Important Google Drive Documents:

Miscellanious commands:
https://docs.google.com/document/d/1YIPwidYIA1nwAA2GbczWpTUEh0doGIe9/edit?usp=sharing&ouid=118151059404127538464&rtpof=true&sd=true

# Install MATLAB:
Install statistics and machine learning toolbox

# Other Github Repos Used in this Project:

### https://github.com/ncbi/sra-tools
Description: Used for downloading transcriptomes in fastq format from the SRA (sequence read archive).

### https://github.com/alexdobin/STAR/tree/master
Description: Use for aligning RNA-seq datasets to HiC assembly. See miscellanious commands document.

# Descriptions of code:

### colour_TPM_using_bed.js
Description: Used to load divergent colourmap for visualizing strand-specific transcription levels in 3D space using CSynth GUI.

### chromosome_3D_alignment.m
Description: Used to align cell cycle stage specific structures from the same chromosome using a rigid rotation.

### group_TADs.m
Description: Used to group dinoflagellate Topological Associating Domains (TADs) into groups based on their CSynth asphericity.

### increase_contact_strength.m
Description: Used to control for differences in unligated (Unligated DNA pieces look like two contacts with a separation of a few hundred bp) and self-ligated (formed when two ends of the same molecule get ligated and then the formed circle gets cleaved elsewhere) Hi-C artefactual contacts. Unligated and self-ligated molecules do not contain information about the 3D organization of the chromosome. Unligated and Self-ligated molecules are short because digestion resolution <<< HiC resolution. These short molecules can map entirely within a bin, or across a bin boundary. These artefacts affect the value of the main diagonal and 2nd diagonal, and the values of entire matrix after normalization. To ensure off diagonal contacts after normalization are not very small (low force & slow convergence in CSynth molecular dynamics run), lower value of main diagonal and 2nd diagonal prior to Hi-C matrix balancing. This ensures the CSynth structure computation converges rapidly.

### sim_HiC_map_CSynth.m
Description: Used to simulate Hi-C contact map from CSynth structure. Compare this to original Hi-C map for an assessment of CSynth accuracy.

### Ps_from_CSynth.m
Description: Used to calculate contact probability "P(s)" curves from CSynth structures.

### percentage_chr_unmappable_vs_read-length.m
Description: Estimate the % of unmappable sequence in the chromosome using reads of variable length by calculating the frequency of strings of consecutive ambiguous bases of different length.

### RNAseq_align_to_CSynth_agg_chr.m
Description: Spatial patterns of active transcription determined by mapping RNA-seq data to 3D position on CSynth structure. An aggregate analysis is done where chromosomes are isotropically scaled to the same size, and integrated in cylindrical coordinates to show the average transcription level vs distance to centre/surface of the chromosome. Consider adding 1D line plot.

### Spatial_Nematic_Order.m
Description: Used to calculate alignment (nematic order parameter) of tangent vectors to DNA that are in close spatial proximity.

### TAD_analysis.m
Description: Used to calculate average level and strand of transcription near TAD boundary.

### active_gene_density_corr_TAD_asphericity.m
Description: Calculate density of active genes within TAD and correlate this with 3D TAD asphericity as determined by CSynth structure.

### convert_hic_to_CSynth.py
Description: A way to directly convert .hic formatted Hi-C data to a format accepted by CSynth.

### chrom_length_order.py
Description: Finds the largest scaffolds or "chromosomes" from dinoflagellate HiC data. Suprisingly, scaffolds are not ordered perfectly by size.

### CSynth_FISH_Compare_Su_Cell_2020.m
Description: Use to assess accuracy of CSynth conformations and to optomize CSynth parameters, using a human cell line where HiC data and 3D FISH data (652 probes) exists.

### CSynth_FISH_Compare_Wang_Science_2016.m
Description: Older version of CSynth_FISH_Compare_Su_Cell_2020.m using a different FISH data set with fewer probes (30).

### dinoflagellate_Ps.m
Description: Possibly useless. Used to compute contact probability "P(s)" curves from dinoflagellate HiC data. Use the python cooltools package instead.

### fractal_equilbrium_load.m
Description: Used to compute contact probability "P(s)" curves from fractal and equilibrium globule conformations as a positive control for the calculation.

### Cholesteric_HiC.m
Description: Generate cholesteric polymer model with extra chromosomal loops. Vary loop length, number of discs, cholesteric pitch, and other model parameters. Compute cholesteric HiC matrix and contact probability "P(s)" curve.

### iter_Cholesteric_HiC.m
Description: Generate cholesteric polymer model with extra chromosomal loops. Vary loop length, number of discs, cholesteric pitch, and other model parameters. Compute cholesteric HiC matrix and contact probability "P(s)" curve.

### constrained_RW_1D.m
Description: Function used to generate extrachromosomal loops. Random walks are generated separately for each of the spatial three dimensions.

### constrained_self_avoiding_RW_3D.m
Description: Function used to generate extrachromosomal loops. Combines random walks in each dimension to compute a 3D random walk. 

### contact_probability_xyz.m
Description: Function used to compute contact probability "P(s)" curves from a 3D structure. Calculation is different from computing contact probability curve from HiC matrix.

### Orientation_Order_Parameter.py
Description: Used to compute the correlation between two tangent vectors to the 3D chromosome structure separated by a certain primary sequence length, averaged over the entire chromosome.

### plot_P(s)_tree_of_life.py
Description: Used to compute the contact probability "P(s)" curves for many eukaryotic organisms across the tree of life.

### plot_P(s)_dinoflagellate.py
Description: Used to compute the contact probability "P(s)" curves for dinoflagellate HiC data testing affect of chromosome ordering (chrom_length_order.py) and differences in aseembly (Smic1.0 vs Smic1.1N).

# Data Sources:
### IMR90 Cells HiC data:
Rao, S. et al. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell 159, 1665–1680 (2014).

### IMR90 Cells FISH dataset:
Su, J. H., et al. (2020). Genome-scale imaging of the 3D organization and transcriptional activity of chromatin. Cell, 182(6), 1641-1659.

### Fractal Globule & Equilibrium Globule Structures:
Lieberman-Aiden E., et al. Comprehensive mapping of long-range interactions reveals folding principles of the human genome. Science 2009 Oct 9;326(5950):289-93. PMID: 19815776. See: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199.
GEO Accession #: GSE18199

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102740
From Gibcus, J. H., Samejima, K., Goloborodko, A., Samejima, I., Naumova, N., Nuebler, J., ... & Dekker, J. (2018). A pathway for mitotic chromosome formation. Science, 359(6376), eaao6135. \
human (HeLa) cells, not synchronized: \
GSM2745898_17MAY11_HSC02RPAB-PE50_CCHiC-HeLa-NS-R1_23-05-2011_kittlere.5.10000.cool \
human (HeLa) cells, mitosis: \
GSM2745897_18FEB15_PE50_C66B1AC-A_Sample_HiC1-aka-CCHiC-HeLaS3CCL2p2-M-98.1000.cool

From: Marinov, G. K., Trevino, A. E., Xiang, T., Kundaje, A., Grossman, A. R., & Greenleaf, W. J. (2021). Transcription-dependent domain-scale three-dimensional genome organization in the dinoflagellate Breviolum minutum. Nature Genetics, 53(5), 613-617. \
GSM4658994	Dinoflagellate Breviolum minutum

From: Nand, A., Zhan, Y., Salazar, O. R., Aranda, M., Voolstra, C. R., & Dekker, J. (2021). Genetic and spatial organization of the unusual chromosomes of the dinoflagellate Symbiodinium microadriaticum. Nature Genetics, 53(5), 618-629. \
GSE152150	Dinoflagellate Symbiodinium microadriaticum combined Smic 1.1N


### HiC Data Symbiodinium kawagutii:
HiC Raw Reads
SRR25948349
SRR25948348

Genome Assembly:
http://sampgr.org.cn/index.php/download (V3)