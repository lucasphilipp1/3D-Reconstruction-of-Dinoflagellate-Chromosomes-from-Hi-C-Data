# This GitHub repository contains the code to recreate the analyses of:
Philipp, L., Marinov G. K., Todd S., Weber S. C., 3D Reconstruction of Dinoflagellate Chromosomes from Hi-C Data Refutes the Cholesteric Liquid Crystal Hypothesis. In review.
![auto_symbiodinium_microadriaticum_chr1_3D xyz](https://github.com/user-attachments/assets/1d8cd915-b809-495d-b877-2e65a52e2fd5)
bioRxiv link to pre-print: [https://www.biorxiv.org/content/10.1101/2025.01.24.634729v1](https://www.biorxiv.org/content/10.1101/2025.01.24.634729v1)

# Data Availability:
S. kawagutii Hi-C data, S. kawagutii & S. microadriaticum CSynth structures, aligned RNA-seq .bed files, and TAD .bed files, have been deposited to: https://doi.org/10.5281/zenodo.14285613.

# Descriptions of code:

### CSynthSerial.js
Description: save .png images for a set of CSynth structures. Drag & drop this file onto CSynth first, then drag & drop a set of .xyz files. Images are saved to downloads folder.

### CSynth_FISH_Compare_Su_Cell_2020.m
Description: Use to assess accuracy of CSynth conformations and to optomize CSynth parameters, using a human cell line where Hi-C data and 3D FISH data (652 probes) exist. See: Fig S3 C,D,E,F in paper.

### Fig5_A_CLC_Expression.m
Description: code to create 3D visual model of surface-localized gene expression on CLC chromosomes using a divergent strand-specific expression colormap. See: Fig 5A in paper.

### Orientation_Order_Parameter.py
Description: Used to compute the average correlation between two tangent vectors to the 3D chromosome structure separated by a given amount of primary sequence, averaged over the entire chromosome. See: Fig 4B in paper. Code adapted from: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.248101.

### RNAseq_align_to_CSynth_agg_chr.m
Description: Visualization of spatial variation in transcription levels determined by mapping RNA-seq data to 3D position on CSynth structure. An aggregate analysis is done where chromosomes are isotropically scaled to the same size, viewed using cylindrical coordinates to show the average transcription level vs distance to centre/surface of the chromosome. See: Fig 5D & S9 in paper.

### Spatial_Nematic_Order.m
Description: Used to calculate alignment (nematic order parameter) of tangent vectors to DNA that are in close spatial proximity. See: Fig 3C in paper.
![nematic_order_param](https://github.com/user-attachments/assets/bb09b6e8-f223-44b6-b24a-973050fcd4e6)

### TAD_analysis.m
Description: Used to calculate average level and strand of transcription near TAD boundary. See: Fig 3C in paper for histograms of TAD asphericity. See: Fig S10 in paper for convergent transcription at TAD boundaries.

### chr_coverage_at_various_TPM.m
Description: Cumulative chromosomal coverage of active gene sequence with TPM > x. This coverage is compared to the chromosome surface area to volume ratio. See Fig S8 in paper.

### colour_TPM_using_bed.js
Description: Used to load divergent colourmap for visualizing strand-specific transcription levels in 3D space using CSynth GUI. Drag and drop this file into CSynth browswer window before or after loading Hi-C contact data, as in Fig 3B.

![RNAseq_colormap](https://github.com/user-attachments/assets/a1b036c4-f829-4ea9-969f-de1760b14e9d)
[<img align="right" width="100" height="100" src="<img width="880" alt="RNAseq+CSynth" src="">](https://github.com/user-attachments/assets/2b0820df-3fa7-4f30-9a2a-6024b9f9bfc6)

### constrained_RW_1D.m
Description: Function used to generate extrachromosomal loops. Random walks are generated separately for each of the spatial three dimensions.

### constrained_self_avoiding_RW_3D.m
Description: Function used to generate extrachromosomal loops. Combines random walks in each dimension to compute a 3D random walk.

### contact_probability_xyz.m
Description: Function used to compute contact probability curves from a 3D structure. Calculation is different from computing contact probability curve from Hi-C matrix. See: Naumova, N., Imakaev, M., Fudenberg, G., Zhan, Y., Lajoie, B. R., Mirny, L. A., & Dekker, J. (2013). Organization of the mitotic chromosome. Science, 342(6161), 948-953. for detailed methods.

### convert_hic_to_CSynth.py
Description: Convert .hic formatted Hi-C data to a format readable by CSynth. Create separate files for each chromosome. Column 1: primary sequence location [bp], Column 2: primary sequence location [bp], Column 3: Hi-C contact strength. Missing values are rare and likely due to repetative regions in the assembly.

<img width="183" alt="Screenshot 2024-10-14 at 3 28 36 PM" src="https://github.com/user-attachments/assets/9855140a-cfe3-4e90-b4b7-09638ab92f7f">

### convert_cool_to_CSynth.sh
Description: Convert .cool or .mcool formatted Hi-C data to a format readable by CSynth. Create separate files for each chromosome. Run in a new terminal window after installing the cooler package: https://github.com/open2c/cooler. 

### cool_plot_HiC.m
Description: Visualize Hi-C contact maps using the cooltools package. A minor modification of:
https://cooltools.readthedocs.io/en/latest/notebooks/contacts_vs_distance.html. See: Fig 1B in the paper.

### fractal_equilbrium_load.m
Description: Used to compute contact probability curves from fractal and equilibrium globule as a positive control to validate the code later applied to CLC model chromosomes.

### iter_Cholesteric_HiC.m
Description: Generate model cholesteric liquid crystal (CLC) chromosomes with extra chromosomal loops. Many adjustable parameters including: loop length, number of discs, cholesteric pitch, and more. Compute single cell and population-level HiC matrices and contact probability curves. See Fig. 1A, 2A, S1, S2, S6 in the paper.

<img width="1118" alt="figure_S6" src="https://github.com/user-attachments/assets/ff609024-7e67-457a-a055-ea41388e2521" />

### plot_P(s)_dinoflagellate.py
Description: Used to compute contact probability curves from dinoflagellate Hi-C data. See: Fig 2B, C, & D in the paper.

### sim_HiC_map_CSynth.m
Description: Used to simulate IMR90 chr21 Hi-C contact map from CSynth structure. This contact map is compared with empirical Hi-C contact for an assessment of CSynth accuracy. See: Fig S3 A & G in paper.

# Data Sources:
### Hi-C data pre-processing:
For an overview of how dinoflagellate contact maps are assembled from raw reads see (Marinov, G., et al. 2024). Previously published Hi-C assisted genome assemblies and contact strengths were accessed at: GEO accession: GSE152150, GSE152150_HiC-Dplus.smic1.1N.mapq_30.1000.mcool (Nand, A., et al., 2021), https://doi.org/10.5281/zenodo.10035644 (Marinov, G., et al. 2024), and GSE153950 (Marinov, G., et al. 2021).

### Human control Hi-C and FISH data:
Previously published human IMR90 Hi-C data (Rao, S. S., et al. 2014) was accessed from the GEO using accession number: GSE63525. Previously published IMR90 FISH probe 3D locations (Su, J. H., et al. 2020) were downloaded from: https://doi.org/10.5281/zenodo.3928890.

### Equilibrium globule structures:
Equilibrium globule structures were downloaded from the GEO using the accession number: GSE18199 (Lieberman-Aiden, E., et al., 2009).

### Spatial Organization of Transcription Analysis:
RNA-seq data used in this study can be accessed via the SRA using the accession numbers:
SRR3337493 (S. microadriaticum, (Liew, Y. J., et al., 2017)), SRR9417753 SRR9417755 SRR9417756 SRR1300302 SRR1300303 SRR1300304 SRR1300305 (S. kawagutii, (Li, T., et al., 2020; Keeling, P. J., et al., 2014)).

# Questions:
If you have questions about this repository please contact Lucas Philipp (lucas.philipp@mail.mcgill.ca).
