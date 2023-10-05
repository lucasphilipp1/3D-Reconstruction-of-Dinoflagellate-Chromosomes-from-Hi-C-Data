clc
clear
chromosome = importdata('symbiodinium_microadriaticum_chr2.xyz');
RNAseq = importdata('s_microadriaticum_chr2_pilon_strandsum_TPM.bed');

HiC_resolution = chromosome(2,1)-chromosome(1,1);

genes_start_end_TPM = [RNAseq.data(:,1) RNAseq.data(:,2) RNAseq.data(:,8)];

[coeff,score] = pca(chromosome(:,2:4));

chromosome_PCA = [chromosome(:,1) score];

%get xyz coordinates for every base pair by linearly interpolating between the model
chromosome_PCA_interpolated = zeros(chromosome(end,1),4);
chromosome_PCA_interpolated(:,1) = linspace(1,chromosome(end,1),chromosome(end,1))';
for i = 1:1:size(chromosome_PCA,1)-1
        chromosome_PCA_interpolated(HiC_resolution*(i-1)+1:HiC_resolution*i,2:4) = [linspace(chromosome_PCA(i,2),chromosome_PCA(i+1,2),HiC_resolution)' linspace(chromosome_PCA(i,3),chromosome_PCA(i+1,3),HiC_resolution)' linspace(chromosome_PCA(i,4),chromosome_PCA(i+1,4),HiC_resolution)'];
end

genes_xyz = [];
genes_center_xyz = [];

%store coordinates for each base pair in each gene
for i = 1:1:size(genes_start_end_TPM,1)
    genes_xyz = [genes_xyz; chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),2) chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),3) chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),4)];
    genes_center_xyz = [genes_center_xyz; mean(chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),2)) mean(chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),3)) mean(chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),4))];
    size(genes_start_end_TPM,1)-i
end

%column 1 is PC1, the chromosome long axis
[theta_TPM,rho_TPM z_TPM] = cart2pol(genes_center_xyz(:,2),genes_center_xyz(:,3),genes_center_xyz(:,1));

figure
scatter(z_TPM,rho_TPM,[],RNAseq.data(:,8), 'filled')
xlim([-100 100])
ylim([0 100])
xlabel('cylidrical axis','FontSize', 24)
ylabel('radial axis','FontSize', 24)
cb = colorbar; 
caxis([min(RNAseq.data(:,8)),max(RNAseq.data(:,8))]);
cb.Label.String = 'TPM';
cb.FontSize = 16;

%gene density is enriched 2.5Mbs away from the telomere, measured in 100-kb windows

%genes per unit volume plot
chrom_xyz = chromosome_PCA_interpolated(:,2:4);

%normalized gene density
%% Compute density
% Put points into 3D bins; xyzBinNum is an nx3 matrix containing
% the bin ID for n values in xyz for the [x,y,z] axes.
nBins = 10;  % number of bins
chrom_xbins = linspace(min(chrom_xyz(:,1)),max(chrom_xyz(:,1))*1,nBins+1);
chrom_ybins = linspace(min(chrom_xyz(:,2)),max(chrom_xyz(:,2))*1,nBins+1);
chrom_zbins = linspace(min(chrom_xyz(:,3)),max(chrom_xyz(:,3))*1,nBins+1);
chrom_xyzBinNum = [...
    discretize(chrom_xyz(:,1),chrom_xbins), ...
    discretize(chrom_xyz(:,2),chrom_ybins), ...
    discretize(chrom_xyz(:,3),chrom_zbins), ...
    ];

genes_xyzBinNum = [...
    discretize(genes_xyz(:,1),chrom_xbins), ...
    discretize(genes_xyz(:,2),chrom_ybins), ...
    discretize(genes_xyz(:,3),chrom_zbins), ...
    ];
% bin3D is a mx3 matrix of m unique 3D bins that appear 
% in xyzBinNum, sorted.  binNum is a nx1 vector of bin
% numbers identifying the bin for each xyz point. For example,
% b=xyz(j,:) belongs to bins3D(b,:).
[chrom_bins3D, ~, chrom_binNum] = unique(chrom_xyzBinNum, 'rows');

%map gene_density onto the same bins as chrom_density
gene_density = [];
for i = 1:1:size(chrom_bins3D,1)
        gene_density = [gene_density; sum(ismember(genes_xyzBinNum, chrom_bins3D(i,:),'rows'));];
        size(chrom_bins3D,1)-i
end

% density is a mx1 vector of integers showing the number of 
% xyz points in each of the bins3D. To see the number of points
% in bins3D(k,:), density(k).  
chrom_density = histcounts(chrom_binNum,[1:size(chrom_bins3D,1),inf])'; 

% Compute bin centers
chrom_xbinCnt = chrom_xbins(2:end)-diff(chrom_xbins)/2;
chrom_ybinCnt = chrom_ybins(2:end)-diff(chrom_ybins)/2;
chrom_zbinCnt = chrom_zbins(2:end)-diff(chrom_zbins)/2;

chrom_bins3D(find(gene_density==0),:) = [];
chrom_density(find(gene_density==0)) = [];
gene_density(find(gene_density==0)) = [];

%% Plot raw data
fig = figure();
%% Plot scatter3
scatter3(...
    chrom_xbinCnt(chrom_bins3D(:,1)), ...
    chrom_ybinCnt(chrom_bins3D(:,2)), ...
    chrom_zbinCnt(chrom_bins3D(:,3)), ...
    gene_density./chrom_density*100, ...
    gene_density./chrom_density, 'filled', ...
    'MarkerFaceAlpha',.6)
axis equal
box on
xlabel('x')
ylabel('y')
zlabel('z')
cb = colorbar
cb.Label.String = 'Gene Density/Chromsome Density';
cb.FontSize = 16;
