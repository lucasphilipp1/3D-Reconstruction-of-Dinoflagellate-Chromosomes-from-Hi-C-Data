# This GitHub repository contains the code to recreate the analyses of:
Philipp, L., Marinov G., Todd S., Weber S., 3D Reconstruction of Dinoflagellate Chromosomes from Hi-C Data Refutes the Cholesteric Liquid Crystal Hypothesis. In review.

bioRxiv link to pre-print: __________

# Zenodo
Symbiodinium kawagutii and Symbiodinium microadriaticum chromosome structures have been deposited at: 

# Descriptions of code:

### colour_TPM_using_bed.js
Description: Used to load divergent colourmap for visualizing strand-specific transcription levels in 3D space using CSynth GUI.

### sim_HiC_map_CSynth.m
Description: Used to simulate Hi-C contact map from CSynth structure. Compare this to original Hi-C map for an assessment of CSynth accuracy.

%Fig S3 A & G in paper

### RNAseq_align_to_CSynth_agg_chr.m
Description: Spatial patterns of active transcription determined by mapping RNA-seq data to 3D position on CSynth structure. An aggregate analysis is done where chromosomes are isotropically scaled to the same size, and integrated in cylindrical coordinates to show the average transcription level vs distance to centre/surface of the chromosome. Consider adding 1D line plot.

%Fig 3D & S9 in paper
%visualization of spatial variation in transcription levels using a
%cylindrical coorindate system

### Spatial_Nematic_Order.m
Description: Used to calculate alignment (nematic order parameter) of tangent vectors to DNA that are in close spatial proximity.
%Fig 3 C in paper

### TAD_analysis.m
Description: Used to calculate average level and strand of transcription near TAD boundary.

%Fig 4 B in paper: histograms of TAD asphericity
%Fig S11 in paper: convergent transcription at TAD boundaries

### TAD_asphericity_correlations.m
Description:

%Fig 4 D,E,F in paper:
%D) Active gene density is not correlated with TAD asphericity. 
%E) Transcription magnitude (transcripts per million) is not correlated with TAD asphericity.
%F) The length of active gene arrays is not correlated with TAD asphericity.

%Fig S12 B, histograph of TAD sizes

### convert_hic_to_CSynth.py
Description: A way to directly convert .hic formatted Hi-C data to a format accepted by CSynth.

### CSynth_FISH_Compare_Su_Cell_2020.m
Description: Use to assess accuracy of CSynth conformations and to optomize CSynth parameters, using a human cell line where HiC data and 3D FISH data (652 probes) exists.
% Fig S3 C,D,E,F in paper

### fractal_equilbrium_load.m
Description: Used to compute contact probability "P(s)" curves from fractal and equilibrium globule conformations as a positive control for the calculation.

### iter_Cholesteric_HiC.m
Description: Generate cholesteric polymer model with extra chromosomal loops. Vary loop length, number of discs, cholesteric pitch, and other model parameters. Compute cholesteric HiC matrix and contact probability curves.

### constrained_RW_1D.m
Description: Function used to generate extrachromosomal loops. Random walks are generated separately for each of the spatial three dimensions.

### constrained_self_avoiding_RW_3D.m
Description: Function used to generate extrachromosomal loops. Combines random walks in each dimension to compute a 3D random walk. 

### contact_probability_xyz.m
Description: Function used to compute contact probability "P(s)" curves from a 3D structure. Calculation is different from computing contact probability curve from HiC matrix.

### Orientation_Order_Parameter.py
Description: Used to compute the correlation between two tangent vectors to the 3D chromosome structure separated by a certain primary sequence length, averaged over the entire chromosome.

### plot_P(s)_dinoflagellate.py
Description: Used to compute the contact probability "P(s)" curves for dinoflagellate HiC data testing affect of chromosome ordering (chrom_length_order.py) and differences in aseembly (Smic1.0 vs Smic1.1N).

### Fig3_A_CLC_Expression.m
Description: code to create 3D visual model of surface-localized gene expression on CLC chromosomes using a divergent strand-specific expression colormap

### chr_coverage_at_various_TPM.m
%Fig S10
%Cumulative chromosomal coverage of active gene sequence with TPM > x. 
%chromosome surface area to volume ratio

### asphericity_vs_size.m
%Fig S13 in the paper
%quality control of asphericity metrics
%calculate asphericity for randomly selected regions of various sizes

# Data Sources:
### IMR90 Cells HiC data:
Rao, S. et al. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell 159, 1665–1680 (2014).

### IMR90 Cells FISH dataset:
Su, J. H., et al. (2020). Genome-scale imaging of the 3D organization and transcriptional activity of chromatin. Cell, 182(6), 1641-1659.

### Fractal Globule & Equilibrium Globule Structures:
Lieberman-Aiden E., et al. Comprehensive mapping of long-range interactions reveals folding principles of the human genome. Science 2009 Oct 9;326(5950):289-93. PMID: 19815776. See: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199.
GEO Accession #: GSE18199

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

Hi-C data pre-processing:
For an overview of how dinoflagellate contact maps are assembled from raw reads see (Marinov, G., et al. 2024). Previously published Hi-C assisted genome assemblies and contact strengths were accessed at: GEO accession: GSE152150, GSE152150_HiC-Dplus.smic1.1N.mapq_30.1000.mcool (Nand, A., et al., 2021), https://doi.org/10.5281/zenodo.10035644 (Marinov, G., et al. 2024) (check to see that link is working), and GSE153950 (Marinov, G., et al. 2021).
![image](https://github.com/user-attachments/assets/a1750e38-96a1-44f8-a8b6-c16ec75fafd7)

Human control Hi-C and FISH data:
Previously published human IMR90 Hi-C data (Rao, S. S., et al. 2014) was accessed from the GEO using accession number: GSE63525. Previously published IMR90 FISH probe 3D locations (Su, J. H., et al. 2020) were downloaded from: https://doi.org/10.5281/zenodo.3928890.

Equilibrium globule structures:
Equilibrium globule structures were downloaded from the GEO using the accession number: GSE18199 (Lieberman-Aiden, E., et al., 2009).

Spatial Organization of Transcription Analysis:
RNA-seq data used in this study can be accessed via the SRA using the accession numbers:
SRR3337493 (S. microadriaticum, (Liew, Y. J., et al., 2017)), SRR9417753 SRR9417755 SRR9417756 SRR1300302 SRR1300303 SRR1300304 SRR1300305 (S. kawagutii, (Li, T., et al., 2020; Keeling, P. J., et al., 2014)).

Alignment to the genome assembly was performed following the methods outlined in (Marinov, G., et al., 2021).
![image](https://github.com/user-attachments/assets/6e8fbd3b-3ed2-474f-87f8-ffa4e78209b3)


# Questions
If you have questions about this repository please contact Lucas Philipp (lucas.philipp@mail.mcgill.ca).
