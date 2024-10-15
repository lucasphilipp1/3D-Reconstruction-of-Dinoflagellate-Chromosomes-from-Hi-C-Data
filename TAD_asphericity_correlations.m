%Fig 4 D,E,F in paper:
%D) Active gene density is not correlated with TAD asphericity. 
%E) Transcription magnitude (transcripts per million) is not correlated with TAD asphericity.
%F) The length of active gene arrays is not correlated with TAD asphericity.

%Fig S12 B, histograph of TAD sizes

clc
clear

num_chroms = 50;

%%% OLS regression: gene density vs asphericity
TADs_microadriaticum_asphericity_all = importdata('symbiodinium_microadriaticum_allchroms_TADs_asphericity_tensor.txt');
%TADs_microadriaticum_asphericity_all = importdata('symbiodinium_microadriaticum_allchroms_TADs_asphericity_PCA.txt');
microadriaticum_chrom_sizes=importdata('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/s_microadriaticum_chrom.sizes');
data = [];

for i=1:num_chroms

    idx = find(contains(TADs_microadriaticum_asphericity_all.textdata,sprintf('chr%i_pilon',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_microadriaticum_asphericity_all.data(idx(j),1);
            TAD_end = TADs_microadriaticum_asphericity_all.data(idx(j),2);
            asphericity = TADs_microadriaticum_asphericity_all.data(idx(j),3);

            idx_RNA = find((TAD_start < RNAseq.data(:,1) & RNAseq.data(:,1) < TAD_end) | (TAD_start < RNAseq.data(:,2) & RNAseq.data(:,2) < TAD_end)); %find genes in or nearby TADs

            if size(idx_RNA,1)>0

                average_TPM = mean(abs(RNAseq.data(idx_RNA,8)));

                gene_array_edges = [];
                gene_array_length = [];

                temp_RNAseq = RNAseq.data(idx_RNA,:);
                gene_array_edges = find(diff(sign(temp_RNAseq(:,8))));
                gene_array_edges = [1; gene_array_edges; size(temp_RNAseq,1)];

                for k=1:1:size(gene_array_edges,1)-1
                    gene_array_length = [gene_array_length; sum(temp_RNAseq(gene_array_edges(k+1))-temp_RNAseq(gene_array_edges(k)))];
                end
                avg_gene_array_length = mean(abs(gene_array_length));
            else
                average_TPM = 0;
                avg_gene_array_length = 0;
            end
            data = [data; asphericity sum(RNAseq.data(idx_RNA,2)-RNAseq.data(idx_RNA,1)) TAD_end-TAD_start average_TPM avg_gene_array_length];
        end
    else
        % File does not exist.
    end
end

idx_outliers = find(data(:,1)>0.95); %remove outliers
data(idx_outliers,:) = [];

figure
scatter(data(:,2)./data(:,3),data(:,1), 50,[0 0.4470 0.7410],"filled")
ylim([0 0.4])
xlim([0 1])
ax = gca;
ax.FontSize = 16;
ylabel('TAD Asphericity (tensor)', 'fontsize', 24)
xlabel('(Sum Length of Active Genes/TAD Size) [bp/bp]', 'fontsize', 24)

figure
scatter(data(:,3),data(:,1), 50,[0 0.4470 0.7410],"filled")
%ylim([0 1])
ylim([0 0.4])
ax = gca;
ax.FontSize = 16;
ylabel('TAD Asphericity (tensor)', 'fontsize', 24)
xlabel('TAD Size [bp]', 'fontsize', 24)

figure
scatter(data(:,4),data(:,1), 50,[0 0.4470 0.7410],"filled")
%ylim([0 1])
ylim([0 0.4])
xlim([0 50])
ax = gca;
ax.FontSize = 16;
ylabel('TAD Asphericity (tensor)', 'fontsize', 24)
xlabel('Average TPM in TAD', 'fontsize', 24)

figure
scatter(data(:,5)./data(:,3),data(:,1), 50,[0 0.4470 0.7410],"filled")
%ylim([0 1])
ylim([0 0.4])
xlim([0 1])
ax = gca;
ax.FontSize = 16;
ylabel('TAD Asphericity (tensor)', 'fontsize', 24)
xlabel('Average Gene Array Length/TAD Size [bp/bp]', 'fontsize', 24)

%Fig S12 B, histograph of TAD sizes
figure
histogram(data(:,3),'FaceColor',[0 0.4470 0.7410])
ax = gca;
ax.FontSize = 16;
ylabel('Count', 'fontsize', 24)
xlabel('TAD Size [bp]', 'fontsize', 24)
xlim([0 5*10^6])
ylim([0 80])


TADs_kawagutii_asphericity_all = importdata('symbiodinium_kawagutii_allchroms_TADs_asphericity_tensor.txt');
%TADs_kawagutii_asphericity_all = importdata('symbiodinium_kawagutii_allchroms_TADs_asphericity_PCA.txt');
kawagutii_chrom_sizes=importdata('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Large_Files/Dinoflagellate HiC Data/s_kawagutii_chrom.sizes');
data = [];

for i=1:num_chroms

    idx = find(contains(TADs_kawagutii_asphericity_all.textdata,sprintf('HiC_scaffold_%i',i)));

    if size(idx,1)>0
        % File exists.
        RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));

        for j = 1:1:size(idx,1)
            TAD_start = TADs_kawagutii_asphericity_all.data(idx(j),1);
            TAD_end = TADs_kawagutii_asphericity_all.data(idx(j),2);
            asphericity = TADs_kawagutii_asphericity_all.data(idx(j),3);

            idx_RNA = find((TAD_start < RNAseq.data(:,1) & RNAseq.data(:,1) < TAD_end) | (TAD_start < RNAseq.data(:,2) & RNAseq.data(:,2) < TAD_end)); %find genes in or nearby TADs

            if size(idx_RNA,1)>0

                average_TPM = mean(abs(RNAseq.data(idx_RNA,8)));

                gene_array_edges = [];
                gene_array_length = [];

                temp_RNAseq = RNAseq.data(idx_RNA,:);
                gene_array_edges = find(diff(sign(temp_RNAseq(:,8))));
                gene_array_edges = [1; gene_array_edges; size(temp_RNAseq,1)];

                for k=1:1:size(gene_array_edges,1)-1
                    gene_array_length = [gene_array_length; sum(temp_RNAseq(gene_array_edges(k+1))-temp_RNAseq(gene_array_edges(k)))];
                end
                avg_gene_array_length = mean(abs(gene_array_length));
            else
                average_TPM = 0;
                avg_gene_array_length = 0;
            end

            data = [data; asphericity sum(RNAseq.data(idx_RNA,2)-RNAseq.data(idx_RNA,1)) TAD_end-TAD_start average_TPM avg_gene_array_length];
        end
    else
        % File does not exist.
    end
end

idx_outliers = find(data(:,1)>0.95); %remove outliers
data(idx_outliers,:) = [];

figure
scatter(data(:,2)./data(:,3),data(:,1), 50,[0.4660 0.6740 0.1880],"filled")
%ylim([0 1])
ylim([0 0.4])
xlim([0 1])
ax = gca;
ax.FontSize = 16;
ylabel('TAD Asphericity (tensor)', 'fontsize', 24)
xlabel('(Sum Length of Active Genes/TAD Size) [bp/bp]', 'fontsize', 24)

figure
scatter(data(:,3),data(:,1), 50,[0.4660 0.6740 0.1880],"filled")
%ylim([0 1])
ylim([0 0.4])
ax = gca;
ax.FontSize = 16;
ylabel('TAD Asphericity (tensor)', 'fontsize', 24)
xlabel('TAD Size [bp]', 'fontsize', 24)

figure
scatter(data(:,4),data(:,1), 50,[0.4660 0.6740 0.1880],"filled")
%ylim([0 1])
ylim([0 0.4])
xlim([0 50])
ax = gca;
ax.FontSize = 16;
ylabel('TAD Asphericity (tensor)', 'fontsize', 24)
xlabel('Average TPM in TAD', 'fontsize', 24)

figure
scatter(data(:,5)./data(:,3),data(:,1), 50,[0.4660 0.6740 0.1880],"filled")
%ylim([0 1])
ylim([0 0.4])
xlim([0 1])
ax = gca;
ax.FontSize = 16;
ylabel('TAD Asphericity (tensor)', 'fontsize', 24)
xlabel('Average Gene Array Length/TAD Size [bp/bp]', 'fontsize', 24)

%Fig S12 B, histograph of TAD sizes
figure
histogram(data(:,3),'FaceColor',[0.4660 0.6740 0.1880])
ax = gca;
ax.FontSize = 16;
ylabel('Count', 'fontsize', 24)
xlabel('TAD Size [bp]', 'fontsize', 24)
xlim([0 5*10^6])
ylim([0 80])