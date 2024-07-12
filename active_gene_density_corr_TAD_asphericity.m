clc
clear

num_chroms = 52;
dist_to_TAD = 0;

num_genes_in_TAD_asphericity_quartile_kawagutii = {};
gene_density_in_TAD_asphericity_quartile_kawagutii = {};
TPM_genes_in_TAD_asphericity_quartile_kawagutii = {};

num_genes_in_TAD = [];
gene_density_in_TAD = [];
TPM_genes_in_TAD = [];

%ASPHERICITY FIRST QUARTILE

TADs_kawagutii_asphericity_first_quartile = importdata('symbiodinium_kawagutii_allchroms_TADs_asphericity_first_quartile.txt');

for i=1:num_chroms

    idx = find(contains(TADs_kawagutii_asphericity_first_quartile.textdata,sprintf('HiC_scaffold_%i',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_kawagutii_asphericity_first_quartile.data(idx(j),1);
            TAD_end = TADs_kawagutii_asphericity_first_quartile.data(idx(j),2);

            %search for genes within +/- dist_to_TAD bp from TAD start/end
            idx_RNA = find(TAD_start - dist_to_TAD < RNAseq.data(:,1) & RNAseq.data(:,2) < TAD_end + dist_to_TAD); %find genes in or nearby TADs

            num_genes_in_TAD = [num_genes_in_TAD; size(idx_RNA,1)];
            gene_density_in_TAD = [gene_density_in_TAD; size(idx_RNA,1)/(TAD_end-TAD_start)];
            TPM_genes_in_TAD = [TPM_genes_in_TAD; abs(RNAseq.data(idx_RNA,8))];
        end
    else
        % File does not exist.
    end
end

idx_outliers = find(TPM_genes_in_TAD>prctile(TPM_genes_in_TAD,99)); %remove outliers
TPM_genes_in_TAD(idx_outliers) = [];

num_genes_in_TAD_asphericity_quartile_kawagutii{1} = num_genes_in_TAD;
gene_density_in_TAD_asphericity_quartile_kawagutii{1} = gene_density_in_TAD;
TPM_genes_in_TAD_asphericity_quartile_kawagutii{1} = TPM_genes_in_TAD;

num_genes_in_TAD = [];
gene_density_in_TAD = [];
TPM_genes_in_TAD = [];

%ASPHERICITY SECOND QUARTILE

TADs_kawagutii_asphericity_second_quartile = importdata('symbiodinium_kawagutii_allchroms_TADs_asphericity_second_quartile.txt');

for i=1:num_chroms

    idx = find(contains(TADs_kawagutii_asphericity_second_quartile.textdata,sprintf('HiC_scaffold_%i',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_kawagutii_asphericity_second_quartile.data(idx(j),1);
            TAD_end = TADs_kawagutii_asphericity_second_quartile.data(idx(j),2);

            %search for genes within +/- dist_to_TAD bp from TAD start/end
            idx_RNA = find(TAD_start - dist_to_TAD < RNAseq.data(:,1) & RNAseq.data(:,2) < TAD_end + dist_to_TAD); %find genes in or nearby TADs

            num_genes_in_TAD = [num_genes_in_TAD; size(idx_RNA,1)];
            gene_density_in_TAD = [gene_density_in_TAD; size(idx_RNA,1)/(TAD_end-TAD_start)];
            TPM_genes_in_TAD = [TPM_genes_in_TAD; abs(RNAseq.data(idx_RNA,8))];
        end
    else
        % File does not exist.
    end
end

idx_outliers = find(TPM_genes_in_TAD>prctile(TPM_genes_in_TAD,99)); %remove outliers
TPM_genes_in_TAD(idx_outliers) = [];

num_genes_in_TAD_asphericity_quartile_kawagutii{2} = num_genes_in_TAD;
gene_density_in_TAD_asphericity_quartile_kawagutii{2} = gene_density_in_TAD;
TPM_genes_in_TAD_asphericity_quartile_kawagutii{2} = TPM_genes_in_TAD;

num_genes_in_TAD = [];
gene_density_in_TAD = [];
TPM_genes_in_TAD = [];

%ASPHERICITY THIRD QUARTILE

TADs_kawagutii_asphericity_third_quartile = importdata('symbiodinium_kawagutii_allchroms_TADs_asphericity_third_quartile.txt');

for i=1:num_chroms

    idx = find(contains(TADs_kawagutii_asphericity_third_quartile.textdata,sprintf('HiC_scaffold_%i',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_kawagutii_asphericity_third_quartile.data(idx(j),1);
            TAD_end = TADs_kawagutii_asphericity_third_quartile.data(idx(j),2);

            %search for genes within +/- dist_to_TAD bp from TAD start/end
            idx_RNA = find(TAD_start - dist_to_TAD < RNAseq.data(:,1) & RNAseq.data(:,2) < TAD_end + dist_to_TAD); %find genes in or nearby TADs

            num_genes_in_TAD = [num_genes_in_TAD; size(idx_RNA,1)];
            gene_density_in_TAD = [gene_density_in_TAD; size(idx_RNA,1)/(TAD_end-TAD_start)];
            TPM_genes_in_TAD = [TPM_genes_in_TAD; abs(RNAseq.data(idx_RNA,8))];
        end
    else
        % File does not exist.
    end
end

idx_outliers = find(TPM_genes_in_TAD>prctile(TPM_genes_in_TAD,99)); %remove outliers
TPM_genes_in_TAD(idx_outliers) = [];

num_genes_in_TAD_asphericity_quartile_kawagutii{3} = num_genes_in_TAD;
gene_density_in_TAD_asphericity_quartile_kawagutii{3} = gene_density_in_TAD;
TPM_genes_in_TAD_asphericity_quartile_kawagutii{3} = TPM_genes_in_TAD;

num_genes_in_TAD = [];
gene_density_in_TAD = [];
TPM_genes_in_TAD = [];

%ASPHERICITY FOURTH QUARTILE

TADs_kawagutii_asphericity_fourth_quartile = importdata('symbiodinium_kawagutii_allchroms_TADs_asphericity_fourth_quartile.txt');

for i=1:num_chroms

    idx = find(contains(TADs_kawagutii_asphericity_fourth_quartile.textdata,sprintf('HiC_scaffold_%i',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_kawagutii_asphericity_fourth_quartile.data(idx(j),1);
            TAD_end = TADs_kawagutii_asphericity_fourth_quartile.data(idx(j),2);

            %search for genes within +/- dist_to_TAD bp from TAD start/end
            idx_RNA = find(TAD_start - dist_to_TAD < RNAseq.data(:,1) & RNAseq.data(:,2) < TAD_end + dist_to_TAD); %find genes in or nearby TADs

            num_genes_in_TAD = [num_genes_in_TAD; size(idx_RNA,1)];
            gene_density_in_TAD = [gene_density_in_TAD; size(idx_RNA,1)/(TAD_end-TAD_start)];
            TPM_genes_in_TAD = [TPM_genes_in_TAD; abs(RNAseq.data(idx_RNA,8))];
        end
    else
        % File does not exist.
    end
end

idx_outliers = find(TPM_genes_in_TAD>prctile(TPM_genes_in_TAD,99)); %remove outliers
TPM_genes_in_TAD(idx_outliers) = [];

num_genes_in_TAD_asphericity_quartile_kawagutii{4} = num_genes_in_TAD;
gene_density_in_TAD_asphericity_quartile_kawagutii{4} = gene_density_in_TAD;
TPM_genes_in_TAD_asphericity_quartile_kawagutii{4} = TPM_genes_in_TAD;

% figure
% violin(num_genes_in_TAD_asphericity_quartile_kawagutii,'facecolor',[0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880]);
% xticks([1 2 3 4])
% xticklabels({'First Quartile','Second Quartile','Third Quartile', 'Fourth Quartile'})
% xtickangle(45)
% ylabel('Number of Active Genes within TAD')
% xlabel('TAD Asphericity')
% ylim([0 750])
% 
% figure
% violin(TPM_genes_in_TAD_asphericity_quartile_kawagutii,'facecolor',[0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880]);
% xticks([1 2 3 4])
% xticklabels({'First Quartile','Second Quartile','Third Quartile', 'Fourth Quartile'})
% xtickangle(45)
% ylabel('TPM of Active Genes within TAD')
% xlabel('TAD Asphericity')
% ylim([0 100])

figure
violin(gene_density_in_TAD_asphericity_quartile_kawagutii,'facecolor',[0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880]);
xticks([1 2 3 4])
xticklabels({'First Quartile','Second Quartile','Third Quartile', 'Fourth Quartile'})
xtickangle(45)
ylabel('(# of Active Genes/TAD size) [1/bp]')
xlabel('TAD Asphericity')
ylim([0 5*10^-4])

%%%

num_genes_in_TAD_asphericity_quartile_microadriaticum = {};
gene_density_in_TAD_asphericity_quartile_microadriaticum = {};
TPM_genes_in_TAD_asphericity_quartile_microadriaticum = {};

num_genes_in_TAD = [];
gene_density_in_TAD = [];
TPM_genes_in_TAD = [];

%ASPHERICITY FIRST QUARTILE

TADs_microadriaticum_asphericity_first_quartile = importdata('symbiodinium_microadriaticum_allchroms_TADs_asphericity_first_quartile.txt');

for i=1:num_chroms

    idx = find(contains(TADs_microadriaticum_asphericity_first_quartile.textdata,sprintf('chr%i_pilon',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_microadriaticum_asphericity_first_quartile.data(idx(j),1);
            TAD_end = TADs_microadriaticum_asphericity_first_quartile.data(idx(j),2);

            %search for genes within +/- dist_to_TAD bp from TAD start/end
            idx_RNA = find(TAD_start - dist_to_TAD < RNAseq.data(:,1) & RNAseq.data(:,2) < TAD_end + dist_to_TAD); %find genes in or nearby TADs

            num_genes_in_TAD = [num_genes_in_TAD; size(idx_RNA,1)];
            gene_density_in_TAD = [gene_density_in_TAD; size(idx_RNA,1)/(TAD_end-TAD_start)];
            TPM_genes_in_TAD = [TPM_genes_in_TAD; abs(RNAseq.data(idx_RNA,8))];
        end
    else
        % File does not exist.
    end
end

idx_outliers = find(TPM_genes_in_TAD>prctile(TPM_genes_in_TAD,99)); %remove outliers
TPM_genes_in_TAD(idx_outliers) = [];

num_genes_in_TAD_asphericity_quartile_microadriaticum{1} = num_genes_in_TAD;
gene_density_in_TAD_asphericity_quartile_microadriaticum{1} = gene_density_in_TAD;
TPM_genes_in_TAD_asphericity_quartile_microadriaticum{1} = TPM_genes_in_TAD;

num_genes_in_TAD = [];
gene_density_in_TAD = [];
TPM_genes_in_TAD = [];

%ASPHERICITY SECOND QUARTILE

TADs_microadriaticum_asphericity_second_quartile = importdata('symbiodinium_microadriaticum_allchroms_TADs_asphericity_second_quartile.txt');

for i=1:num_chroms

    idx = find(contains(TADs_microadriaticum_asphericity_second_quartile.textdata,sprintf('chr%i_pilon',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_microadriaticum_asphericity_second_quartile.data(idx(j),1);
            TAD_end = TADs_microadriaticum_asphericity_second_quartile.data(idx(j),2);

            %search for genes within +/- dist_to_TAD bp from TAD start/end
            idx_RNA = find(TAD_start - dist_to_TAD < RNAseq.data(:,1) & RNAseq.data(:,2) < TAD_end + dist_to_TAD); %find genes in or nearby TADs

            num_genes_in_TAD = [num_genes_in_TAD; size(idx_RNA,1)];
            gene_density_in_TAD = [gene_density_in_TAD; size(idx_RNA,1)/(TAD_end-TAD_start)];
            TPM_genes_in_TAD = [TPM_genes_in_TAD; abs(RNAseq.data(idx_RNA,8))];
        end
    else
        % File does not exist.
    end
end

idx_outliers = find(TPM_genes_in_TAD>prctile(TPM_genes_in_TAD,99)); %remove outliers
TPM_genes_in_TAD(idx_outliers) = [];

num_genes_in_TAD_asphericity_quartile_microadriaticum{2} = num_genes_in_TAD;
gene_density_in_TAD_asphericity_quartile_microadriaticum{2} = gene_density_in_TAD;
TPM_genes_in_TAD_asphericity_quartile_microadriaticum{2} = TPM_genes_in_TAD;

num_genes_in_TAD = [];
gene_density_in_TAD = [];
TPM_genes_in_TAD = [];

%ASPHERICITY THIRD QUARTILE

TADs_microadriaticum_asphericity_third_quartile = importdata('symbiodinium_microadriaticum_allchroms_TADs_asphericity_third_quartile.txt');

for i=1:num_chroms

    idx = find(contains(TADs_microadriaticum_asphericity_third_quartile.textdata,sprintf('chr%i_pilon',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_microadriaticum_asphericity_third_quartile.data(idx(j),1);
            TAD_end = TADs_microadriaticum_asphericity_third_quartile.data(idx(j),2);

            %search for genes within +/- dist_to_TAD bp from TAD start/end
            idx_RNA = find(TAD_start - dist_to_TAD < RNAseq.data(:,1) & RNAseq.data(:,2) < TAD_end + dist_to_TAD); %find genes in or nearby TADs

            num_genes_in_TAD = [num_genes_in_TAD; size(idx_RNA,1)];
            gene_density_in_TAD = [gene_density_in_TAD; size(idx_RNA,1)/(TAD_end-TAD_start)];
            TPM_genes_in_TAD = [TPM_genes_in_TAD; abs(RNAseq.data(idx_RNA,8))];
        end
    else
        % File does not exist.
    end
end

idx_outliers = find(TPM_genes_in_TAD>prctile(TPM_genes_in_TAD,99)); %remove outliers
TPM_genes_in_TAD(idx_outliers) = [];

num_genes_in_TAD_asphericity_quartile_microadriaticum{3} = num_genes_in_TAD;
gene_density_in_TAD_asphericity_quartile_microadriaticum{3} = gene_density_in_TAD;
TPM_genes_in_TAD_asphericity_quartile_microadriaticum{3} = TPM_genes_in_TAD;

num_genes_in_TAD = [];
gene_density_in_TAD = [];
TPM_genes_in_TAD = [];

%ASPHERICITY FOURTH QUARTILE

TADs_microadriaticum_asphericity_fourth_quartile = importdata('symbiodinium_microadriaticum_allchroms_TADs_asphericity_fourth_quartile.txt');

for i=1:num_chroms

    idx = find(contains(TADs_microadriaticum_asphericity_fourth_quartile.textdata,sprintf('chr%i_pilon',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_microadriaticum_asphericity_fourth_quartile.data(idx(j),1);
            TAD_end = TADs_microadriaticum_asphericity_fourth_quartile.data(idx(j),2);

            %search for genes within +/- dist_to_TAD bp from TAD start/end
            idx_RNA = find(TAD_start - dist_to_TAD < RNAseq.data(:,1) & RNAseq.data(:,2) < TAD_end + dist_to_TAD); %find genes in or nearby TADs

            num_genes_in_TAD = [num_genes_in_TAD; size(idx_RNA,1)];
            gene_density_in_TAD = [gene_density_in_TAD; size(idx_RNA,1)/(TAD_end-TAD_start)];
            TPM_genes_in_TAD = [TPM_genes_in_TAD; abs(RNAseq.data(idx_RNA,8))];
        end
    else
        % File does not exist.
    end
end

idx_outliers = find(TPM_genes_in_TAD>prctile(TPM_genes_in_TAD,99)); %remove outliers
TPM_genes_in_TAD(idx_outliers) = [];

num_genes_in_TAD_asphericity_quartile_microadriaticum{4} = num_genes_in_TAD;
gene_density_in_TAD_asphericity_quartile_microadriaticum{4} = gene_density_in_TAD;
TPM_genes_in_TAD_asphericity_quartile_microadriaticum{4} = TPM_genes_in_TAD;

% figure
% violin(num_genes_in_TAD_asphericity_quartile_microadriaticum,'facecolor',[0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410]);
% xticks([1 2 3 4])
% xticklabels({'First Quartile','Second Quartile','Third Quartile', 'Fourth Quartile'})
% xtickangle(45)
% ylabel('Number of Active Genes within TAD')
% xlabel('TAD Asphericity')
% ylim([0 750])
% 
% figure
% violin(TPM_genes_in_TAD_asphericity_quartile_microadriaticum,'facecolor',[0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410]);
% xticks([1 2 3 4])
% xticklabels({'First Quartile','Second Quartile','Third Quartile', 'Fourth Quartile'})
% xtickangle(45)
% ylabel('TPM of Active Genes within TAD')
% xlabel('TAD Asphericity')
% ylim([0 100])

figure
violin(gene_density_in_TAD_asphericity_quartile_microadriaticum,'facecolor',[0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410]);
xticks([1 2 3 4])
xticklabels({'First Quartile','Second Quartile','Third Quartile', 'Fourth Quartile'})
xtickangle(45)
ylabel('(# of Active Genes/TAD size) [1/bp]')
xlabel('TAD Asphericity')
ylim([0 5*10^-4])

%%% OLS regression: gene density vs asphericity 
TADs_microadriaticum_asphericity_all = importdata('symbiodinium_microadriaticum_allchroms_TADs_asphericity.txt');
asphericity_gene_density_in_TAD = [];
for i=1:num_chroms

    idx = find(contains(TADs_microadriaticum_asphericity_all.textdata,sprintf('chr%i_pilon',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_microadriaticum_asphericity_all.data(idx(j),1);
            TAD_end = TADs_microadriaticum_asphericity_all.data(idx(j),2);
            asphericity = TADs_microadriaticum_asphericity_all.data(idx(j),3);

            %search for genes within +/- dist_to_TAD bp from TAD start/end
            idx_RNA = find(TAD_start - dist_to_TAD < RNAseq.data(:,1) & RNAseq.data(:,2) < TAD_end + dist_to_TAD); %find genes in or nearby TADs
            
            asphericity_gene_density_in_TAD = [asphericity_gene_density_in_TAD; asphericity size(idx_RNA,1)/(TAD_end-TAD_start)];

        end
    else
        % File does not exist.
    end
end

figure
scatter(asphericity_gene_density_in_TAD(:,1),asphericity_gene_density_in_TAD(:,2))

b1 = asphericity_gene_density_in_TAD(:,1)\asphericity_gene_density_in_TAD(:,2);

yCalc1 = b1*asphericity_gene_density_in_TAD(:,1);
hold on
plot(asphericity_gene_density_in_TAD(:,1),yCalc1)
xlabel('TAD Asphericity')
ylabel('(# of Active Genes/TAD size) [1/bp]')
title('Symbiodinium microadriatium')

tbl = table(asphericity_gene_density_in_TAD(:,2),asphericity_gene_density_in_TAD(:,1));
mdl = fitlm(tbl,'Var1 ~ Var2');
[p,F] = coefTest(mdl)

TADs_kawagutii_asphericity_all = importdata('s_kawagutii_V3_allchroms_TADs_asphericity.txt');
asphericity_gene_density_in_TAD = [];
for i=1:num_chroms

    idx = find(contains(TADs_kawagutii_asphericity_all.textdata,sprintf('HiC_scaffold_%i',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_kawagutii_asphericity_all.data(idx(j),1);
            TAD_end = TADs_kawagutii_asphericity_all.data(idx(j),2);
            asphericity = TADs_kawagutii_asphericity_all.data(idx(j),3);

            %search for genes within +/- dist_to_TAD bp from TAD start/end
            idx_RNA = find(TAD_start - dist_to_TAD < RNAseq.data(:,1) & RNAseq.data(:,2) < TAD_end + dist_to_TAD); %find genes in or nearby TADs

            asphericity_gene_density_in_TAD = [asphericity_gene_density_in_TAD; asphericity size(idx_RNA,1)/(TAD_end-TAD_start)];

        end
    else
        % File does not exist.
    end
end

idx_outliers = find(asphericity_gene_density_in_TAD(:,1)>0.95); %remove outliers
asphericity_gene_density_in_TAD(idx_outliers,:) = [];

figure
scatter(asphericity_gene_density_in_TAD(:,1),asphericity_gene_density_in_TAD(:,2))

b1 = asphericity_gene_density_in_TAD(:,1)\asphericity_gene_density_in_TAD(:,2);

yCalc1 = b1*asphericity_gene_density_in_TAD(:,1);
hold on
plot(asphericity_gene_density_in_TAD(:,1),yCalc1)
xlabel('TAD Asphericity')
ylabel('(# of Active Genes/TAD size) [1/bp]')
title('Symbiodinium kawagutii')

tbl = table(asphericity_gene_density_in_TAD(:,2),asphericity_gene_density_in_TAD(:,1));
mdl = fitlm(tbl,'Var1 ~ Var2');
[p,F] = coefTest(mdl)


%__________________________________________________________________________
% violin.m - Simple violin plot using matlab default kernel density estimation
% Last update: 10/2015
%__________________________________________________________________________
% This function creates violin plots based on kernel density estimation
% using ksdensity with default settings. Please be careful when comparing pdfs
% estimated with different bandwidth!
%
% Differently to other boxplot functions, you may specify the x-position.
% This is usefule when overlaying with other data / plots.
%__________________________________________________________________________
%
% Please cite this function as:
% Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
% hhoffmann@uni-bonn.de
%
%__________________________________________________________________________
%
% INPUT
%
% Y:     Data to be plotted, being either
%        a) n x m matrix. A 'violin' is plotted for each column m, OR
%        b) 1 x m Cellarry with elements being numerical colums of nx1 length.
%
% varargin:
% xlabel:    xlabel. Set either [] or in the form {'txt1','txt2','txt3',...}
% facecolor: FaceColor. (default [1 0.5 0]); Specify abbrev. or m x 3 matrix (e.g. [1 0 0])
% edgecolor: LineColor. (default 'k'); Specify abbrev. (e.g. 'k' for black); set either [],'' or 'none' if the mean should not be plotted
% facealpha: Alpha value (transparency). default: 0.5
% mc:        Color of the bars indicating the mean. (default 'k'); set either [],'' or 'none' if the mean should not be plotted
% medc:      Color of the bars indicating the median. (default 'r'); set either [],'' or 'none' if the mean should not be plotted
% bw:        Kernel bandwidth. (default []); prescribe if wanted as follows:
%            a) if bw is a single number, bw will be applied to all
%            columns or cells
%            b) if bw is an array of 1xm or mx1, bw(i) will be applied to cell or column (i).
%            c) if bw is empty (default []), the optimal bandwidth for
%            gaussian kernel is used (see Matlab documentation for
%            ksdensity()
%
% OUTPUT
%
% h:     figure handle
% L:     Legend handle
% MX:    Means of groups
% MED:   Medians of groups
% bw:    bandwidth of kernel
%__________________________________________________________________________
%{
% Example1 (default):
disp('this example uses the statistical toolbox')
Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
[h,L,MX,MED]=violin(Y);
ylabel('\Delta [yesno^{-2}]','FontSize',14)
%Example2 (specify facecolor, edgecolor, xlabel):
disp('this example uses the statistical toolbox')
Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
violin(Y,'xlabel',{'a','b','c','d'},'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','b',...
'bw',0.3,...
'mc','k',...
'medc','r--')
ylabel('\Delta [yesno^{-2}]','FontSize',14)
%Example3 (specify x axis location):
disp('this example uses the statistical toolbox')
Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
violin(Y,'x',[-1 .7 3.4 8.8],'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','none',...
'bw',0.3,'mc','k','medc','r-.')
axis([-2 10 -0.5 20])
ylabel('\Delta [yesno^{-2}]','FontSize',14)
%Example4 (Give data as cells with different n):
disp('this example uses the statistical toolbox')
Y{:,1}=rand(10,1);
Y{:,2}=rand(1000,1);
violin(Y,'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','none','bw',0.1,'mc','k','medc','r-.')
ylabel('\Delta [yesno^{-2}]','FontSize',14)
%}
%%
function[h,L,MX,MED,bw]=violin(Y,varargin)
%defaults:
%_____________________
xL=[];
fc=[1 0.5 0];
lc='k';
alp=0.5;
mc='k';
medc='r';
b=[]; %bandwidth
plotlegend=1;
plotmean=1;
plotmedian=1;
x = [];
%_____________________
%convert single columns to cells:
if iscell(Y)==0
    Y = num2cell(Y,1);
end
%get additional input parameters (varargin)
if isempty(find(strcmp(varargin,'xlabel')))==0
    xL = varargin{find(strcmp(varargin,'xlabel'))+1};
end
if isempty(find(strcmp(varargin,'facecolor')))==0
    fc = varargin{find(strcmp(varargin,'facecolor'))+1};
end
if isempty(find(strcmp(varargin,'edgecolor')))==0
    lc = varargin{find(strcmp(varargin,'edgecolor'))+1};
end
if isempty(find(strcmp(varargin,'facealpha')))==0
    alp = varargin{find(strcmp(varargin,'facealpha'))+1};
end
if isempty(find(strcmp(varargin,'mc')))==0
    if isempty(varargin{find(strcmp(varargin,'mc'))+1})==0
        mc = varargin{find(strcmp(varargin,'mc'))+1};
        plotmean = 1;
    else
        plotmean = 0;
    end
end
if isempty(find(strcmp(varargin,'medc')))==0
    if isempty(varargin{find(strcmp(varargin,'medc'))+1})==0
        medc = varargin{find(strcmp(varargin,'medc'))+1};
        plotmedian = 1;
    else
        plotmedian = 0;
    end
end
if isempty(find(strcmp(varargin,'bw')))==0
    b = varargin{find(strcmp(varargin,'bw'))+1}
    if length(b)==1
        disp(['same bandwidth bw = ',num2str(b),' used for all cols'])
        b=repmat(b,size(Y,2),1);
    elseif length(b)~=size(Y,2)
        warning('length(b)~=size(Y,2)')
        error('please provide only one bandwidth or an array of b with same length as columns in the data set')
    end
end
if isempty(find(strcmp(varargin,'plotlegend')))==0
    plotlegend = varargin{find(strcmp(varargin,'plotlegend'))+1};
end
if isempty(find(strcmp(varargin,'x')))==0
    x = varargin{find(strcmp(varargin,'x'))+1};
end
%%
if size(fc,1)==1
    fc=repmat(fc,size(Y,2),1);
end
%% Calculate the kernel density
i=1;
for i=1:size(Y,2)

    if isempty(b)==0
        [f, u, bb]=ksdensity(Y{i},'bandwidth',b(i));
    elseif isempty(b)
        [f, u, bb]=ksdensity(Y{i});
    end

    f=f/max(f)*0.3; %normalize
    F(:,i)=f;
    U(:,i)=u;
    MED(:,i)=nanmedian(Y{i});
    MX(:,i)=nanmean(Y{i});
    bw(:,i)=bb;

end
%%
%-------------------------------------------------------------------------
% Put the figure automatically on a second monitor
% mp = get(0, 'MonitorPositions');
% set(gcf,'Color','w','Position',[mp(end,1)+50 mp(end,2)+50 800 600])
%-------------------------------------------------------------------------
%Check x-value options
if isempty(x)
    x = zeros(size(Y,2));
    setX = 0;
else
    setX = 1;
    if isempty(xL)==0
        disp('_________________________________________________________________')
        warning('Function is not designed for x-axis specification with string label')
        warning('when providing x, xlabel can be set later anyway')
        error('please provide either x or xlabel. not both.')
    end
end
%% Plot the violins
i=1;
for i=i:size(Y,2)
    if isempty(lc) == 1
        if setX == 0
            h(i)=fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
        else
            h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
        end
    else
        if setX == 0
            h(i)=fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
        else
            h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
        end
    end
    hold on
    if setX == 0
        if plotmean == 1
            p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i)) ],[MX(:,i) MX(:,i)],mc,'LineWidth',2);
        end
        if plotmedian == 1
            p(2)=plot([interp1(U(:,i),F(:,i)+i,MED(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i)) ],[MED(:,i) MED(:,i)],medc,'LineWidth',2);
        end
    elseif setX == 1
        if plotmean == 1
            p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i))+x(i)-i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i))+x(i)-i],[MX(:,i) MX(:,i)],mc,'LineWidth',2);
        end
        if plotmedian == 1
            p(2)=plot([interp1(U(:,i),F(:,i)+i,MED(:,i))+x(i)-i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))+x(i)-i],[MED(:,i) MED(:,i)],medc,'LineWidth',2);
        end
    end
end
%% Add legend if requested
if plotlegend==1 & plotmean==1 | plotlegend==1 & plotmedian==1

    if plotmean==1 & plotmedian==1
        L=legend([p(1) p(2)],'Mean','Median');
    elseif plotmean==0 & plotmedian==1
        L=legend([p(2)],'Median');
    elseif plotmean==1 & plotmedian==0
        L=legend([p(1)],'Mean');
    end

    set(L,'box','off','FontSize',14)
else
    L=[];
end
%% Set axis
if setX == 0
    axis([0.5 size(Y,2)+0.5, min(U(:)) max(U(:))]);
elseif setX == 1
    axis([min(x)-0.05*range(x) max(x)+0.05*range(x), min(U(:)) max(U(:))]);
end
%% Set x-labels
xL2={''};
i=1;
for i=1:size(xL,2)
    xL2=[xL2,xL{i},{''}];
end
set(gca,'TickLength',[0 0],'FontSize',12)
box on
if isempty(xL)==0
    set(gca,'XtickLabel',xL2)
end
%-------------------------------------------------------------------------
end %of function
