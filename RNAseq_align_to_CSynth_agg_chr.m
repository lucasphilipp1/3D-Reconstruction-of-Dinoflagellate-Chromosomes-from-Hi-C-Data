clc
clear

num_chroms = 75;
bins = 20;

low_expression_agg = [];
med_expression_agg = [];
high_expression_agg = [];
all_expression_agg = [];

telomeres = [];

P_axis_1 = [];
P_axis_2 = [];
P_axis_3 = [];
asphericity_microadriaticum = [];
asphericity_tensor_microadriaticum =[];

chromosome_length_microadriaticum = [];
mean_TPM_microadriaticum = [];
median_TPM_microadriaticum = [];
percent_gene_content_microadriaticum = [];
mean_intergenic_length_microadriaticum = [];
mean_tandem_gene_array_length_microadriaticum = [];
mean_tandem_gene_array_length_bp_microadriaticum = [];

for i=1:num_chroms
    
    chromosome = importdata(sprintf('symbiodinium_microadriaticum_chr%i_3D.xyz',i));
    RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));

    %chromosome = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i.xyz',i));
    %RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));

    chromosome_length_microadriaticum = [chromosome_length_microadriaticum; chromosome(end,1)];

    chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
    chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));
    chromosome(:,4) = chromosome(:,4) - mean(chromosome(:,4));

    HiC_resolution = chromosome(2,1)-chromosome(1,1);

    genes_start_end_TPM = [RNAseq.data(:,1) RNAseq.data(:,2) RNAseq.data(:,8)];

    mean_TPM_microadriaticum = [mean_TPM_microadriaticum; mean(abs(genes_start_end_TPM(:,3)))];
    median_TPM_microadriaticum = [median_TPM_microadriaticum; median(abs(genes_start_end_TPM(:,3)))];

    percent_gene_content_microadriaticum = [percent_gene_content_microadriaticum; sum(genes_start_end_TPM(:,2)-genes_start_end_TPM(:,1))/chromosome(end,1)*100];
    
    intergenic = [];
    for j = 1:1:size(genes_start_end_TPM,1)
        if j == 1
            intergenic = [intergenic; genes_start_end_TPM(1,1)];
        elseif j == size(genes_start_end_TPM,1)
            intergenic = [intergenic; chromosome(end,1)-genes_start_end_TPM(end,2)];
        else
            intergenic = [intergenic; genes_start_end_TPM(j,1)-genes_start_end_TPM(j-1,2)];
        end
    end
    
    mean_intergenic_length_microadriaticum = [mean_intergenic_length_microadriaticum; mean(intergenic)];
    
    strand_changes = find(diff(sign(genes_start_end_TPM(:,3))))+1;

    tandem_gene_array_lengths = [];
    tandem_gene_array_lengths = [tandem_gene_array_lengths; strand_changes(1)];
    tandem_gene_array_lengths = [tandem_gene_array_lengths; diff(strand_changes)];

    mean_tandem_gene_array_length_microadriaticum = [mean_tandem_gene_array_length_microadriaticum; mean(tandem_gene_array_lengths)];

    strand_changes = [1; strand_changes];

    tandem_gene_array_lengths_bp = [];
    for j = 1:1:size(strand_changes,1)-1
        tandem_gene_array_lengths_bp = [tandem_gene_array_lengths_bp; genes_start_end_TPM(strand_changes(j+1)-1,2)-genes_start_end_TPM(strand_changes(j),1)];
    end

    mean_tandem_gene_array_length_bp_microadriaticum = [mean_tandem_gene_array_length_bp_microadriaticum; mean(tandem_gene_array_lengths_bp)];

    [coeff,score] = pca(chromosome(:,2:4));
    chromosome_PCA = score;

    %calculate moment of inertia tensor
    %Following: Chu, X., & Wang, J. (2023). Quantifying the large-scale chromosome structural dynamics during the mitosis-to-G1 phase transition of cell cycle. Open Biology, 13(11), 230175.
    Itensor=[[chromosome_PCA(:,1)'*chromosome_PCA(:,1) chromosome_PCA(:,1)'*chromosome_PCA(:,2) chromosome_PCA(:,1)'*chromosome_PCA(:,3)]; [chromosome_PCA(:,2)'*chromosome_PCA(:,1) chromosome_PCA(:,2)'*chromosome_PCA(:,2) chromosome_PCA(:,2)'*chromosome_PCA(:,3)]; [chromosome_PCA(:,3)'*chromosome_PCA(:,1) chromosome_PCA(:,3)'*chromosome_PCA(:,2) chromosome_PCA(:,3)'*chromosome_PCA(:,3)]];

    eigenvalues = eig(Itensor);

    asphericity_tensor_microadriaticum = [asphericity_tensor_microadriaticum; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];

    P_axis_1=[P_axis_1; prctile(abs(score(:,1)),95)];
    P_axis_2=[P_axis_2; prctile(abs(score(:,2)),95)];
    P_axis_3=[P_axis_3; prctile(abs(score(:,3)),95)];

    %get xyz coordinates for every base pair by linearly interpolating between the model
    chromosome_PCA_interpolated = zeros(chromosome(end,1),4);
    chromosome_PCA_interpolated(:,1) = linspace(1,chromosome(end,1),chromosome(end,1))';
    for j = 1:1:size(chromosome_PCA,1)-1
        chromosome_PCA_interpolated(HiC_resolution*(j-1)+1:HiC_resolution*j,2:4) = [linspace(chromosome_PCA(j,1),chromosome_PCA(j+1,1),HiC_resolution)' linspace(chromosome_PCA(j,2),chromosome_PCA(j+1,2),HiC_resolution)' linspace(chromosome_PCA(j,3),chromosome_PCA(j+1,3),HiC_resolution)'];
    end

    while genes_start_end_TPM(end,2) > chromosome_PCA_interpolated(end,1)
        genes_start_end_TPM(end,:)=[];
    end

    %genes_xyz =[];
    genes_center_xyz = [];

    %store coordinates for each base pair in each gene
    for k = 1:1:size(genes_start_end_TPM,1)
        %genes_xyz = [genes_xyz; chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),1) chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),2) chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),3)];
        genes_center_xyz = [genes_center_xyz; mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),2)) mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),3)) mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),4))];
    end

    %column 1 is PC1, the chromosome long axis
    %cart2pol does x,y,z -> theta,rho,z. z is unchanged
    [theta_TPM, rho_TPM, z_TPM] = cart2pol(genes_center_xyz(:,2),genes_center_xyz(:,3),genes_center_xyz(:,1));

    %column 1 is PC1, the chromosome long axis
    [theta_chromosome, rho_chromosome, z_chromosome] = cart2pol(chromosome_PCA(:,2),chromosome_PCA(:,3),chromosome_PCA(:,1));

    temp_telomeres = [chromosome_PCA(1,:); chromosome_PCA(end,:)];

    %check that telomeres span a significant amount of clylindrical axis
    %total span of chromosome is [-1 1] (interval of length 2)
    %if  abs(z_chromosome(1,:)./max(abs(z_chromosome)) - z_chromosome(end,:)./max(abs(z_chromosome))) > 1
    if  abs(z_chromosome(1,:)./max(abs(z_chromosome)) - z_chromosome(end,:)./max(abs(z_chromosome))) > 0

        %rescale chromosome dimensions
        all_expression_agg = [all_expression_agg; rho_TPM./max(rho_chromosome) z_TPM./max(abs(z_chromosome)) abs(genes_start_end_TPM(:,3)) ones(size(genes_start_end_TPM,1),1).*i];
        telomeres = [telomeres; rho_chromosome(1,:)./max(rho_chromosome) z_chromosome(1,:)./max(abs(z_chromosome))];
        telomeres = [telomeres; rho_chromosome(end,:)./max(rho_chromosome) z_chromosome(end,:)./max(abs(z_chromosome))];
    end
    i
end

P=prctile(all_expression_agg(:,3),[20 80]);
%low expression <20th percentile
low_expression_agg = all_expression_agg(find(all_expression_agg(:,3)<P(1)),:);

%medium expression
med_expression_agg = all_expression_agg(find(all_expression_agg(:,3)>P(1) & all_expression_agg(:,3)<P(2)),:);

%high expression >80th percentile
high_expression_agg = all_expression_agg(find(all_expression_agg(:,3)>P(2)),:);

figure
[N,c] = hist3(low_expression_agg(:,1:2),[bins bins]);
N=flipud(N);
c{1}=fliplr(c{1});
for i = 1:1:size(N,1)
    N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
end
imagesc(N)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
hold on
%plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k')
hold off
%colormap("sky")
colormap(custom_cmap('green'))
cb = colorbar;
cb.Label.String = 'Normalized Gene Density';
cb.FontSize = 16;
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('cylidrical axis','FontSize', 24)
ylabel('radial axis','FontSize', 24)
t=get(cb,'Limits');
set(cb,'Ticks',[t(1),t(2)])
set(cb,'XTickLabel',{'Low','High',});
%xlim([2.5 20.5])

% figure
% [N,c] = hist3(med_expression_agg(:,1:2),[bins bins]);
% N=flipud(N);
% c{1}=fliplr(c{1});
% for i = 1:1:size(N,1)
%     N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
% end
% imagesc(N)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
% hold on
% plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k')
% hold off
% colormap("sky")
% cb = colorbar;
% cb.Label.String = 'Normalized Gene Density';
% cb.FontSize = 16;
% set(gca,'XTick',[])
% set(gca,'YTick',[])
% xlabel('cylidrical axis','FontSize', 24)
% ylabel('radial axis','FontSize', 24)
% t=get(cb,'Limits');
% set(cb,'Ticks',[t(1),t(2)])
% set(cb,'XTickLabel',{'Low','High',});

figure
[N,c] = hist3(high_expression_agg(:,1:2),[bins bins]);
N=flipud(N);
c{1}=fliplr(c{1});
for i = 1:1:size(N,1)
    N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
end
imagesc(N)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
hold on
%plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k')
hold off
%colormap("sky")
colormap(custom_cmap('green'))
cb = colorbar;
cb.Label.String = 'Normalized Gene Density';
cb.FontSize = 16;
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('cylidrical axis','FontSize', 24)
ylabel('radial axis','FontSize', 24)
t=get(cb,'Limits');
set(cb,'Ticks',[t(1),t(2)])
set(cb,'XTickLabel',{'Low','High',});
%xlim([2.5 20.5])

figure
[N,c] = hist3(all_expression_agg(:,1:2),[bins bins]);
N=flipud(N);
c{1}=fliplr(c{1});
for i = 1:1:size(N,1)
    N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
end
imagesc(N)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
hold on
plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k')
hold off
%colormap("sky")
colormap(custom_cmap('green'))
cb = colorbar;
cb.Label.String = 'Normalized Gene Density';
cb.FontSize = 16;
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('cylidrical axis','FontSize', 24)
ylabel('radial axis','FontSize', 24)
t=get(cb,'Limits');
set(cb,'Ticks',[t(1),t(2)])
set(cb,'XTickLabel',{'Low','High',});
%xlim([2.5 20.5])

ratio_P_axis_1_2_microadriaticum = P_axis_1./P_axis_2;
ratio_P_axis_1_3_microadriaticum = P_axis_1./P_axis_3;
ratio_P_axis_2_3_microadriaticum = P_axis_2./P_axis_3;

for i = 1:1:num_chroms
    asphericity_microadriaticum = [asphericity_microadriaticum; max([ratio_P_axis_1_2_microadriaticum(i), ratio_P_axis_1_3_microadriaticum(i), ratio_P_axis_2_3_microadriaticum(i)])];
end

all_expression_agg_microadriaticum = all_expression_agg;

i=0;
low_expression_agg = [];
med_expression_agg = [];
high_expression_agg = [];
all_expression_agg = [];

telomeres = [];

P_axis_1 = [];
P_axis_2 = [];
P_axis_3 = [];
asphericity_kawagutii = [];
asphericity_tensor_kawagutii=[];

chromosome_length_kawagutii = [];
mean_TPM_kawagutii = [];
median_TPM_kawagutii = [];
percent_gene_content_kawagutii = [];
mean_intergenic_length_kawagutii = [];
mean_tandem_gene_array_length_kawagutii = [];
mean_tandem_gene_array_length_bp_kawagutii = [];

for i=1:num_chroms
    
    %chromosome = importdata(sprintf('symbiodinium_microadriaticum_chr%i_3D.xyz',i));
    %RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));

    chromosome = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i.xyz',i));
    RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));

    chromosome_length_kawagutii = [chromosome_length_kawagutii; chromosome(end,1)];

    chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
    chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));
    chromosome(:,4) = chromosome(:,4) - mean(chromosome(:,4));

    HiC_resolution = chromosome(2,1)-chromosome(1,1);

    genes_start_end_TPM = [RNAseq.data(:,1) RNAseq.data(:,2) RNAseq.data(:,8)];

    mean_TPM_kawagutii = [mean_TPM_kawagutii; mean(abs(genes_start_end_TPM(:,3)))];
    median_TPM_kawagutii = [median_TPM_kawagutii; median(abs(genes_start_end_TPM(:,3)))];

    percent_gene_content_kawagutii = [percent_gene_content_kawagutii; sum(genes_start_end_TPM(:,2)-genes_start_end_TPM(:,1))/chromosome(end,1)*100];
    
    intergenic = [];
    for j = 1:1:size(genes_start_end_TPM,1)
        if j == 1
            intergenic = [intergenic; genes_start_end_TPM(1,1)];
        elseif j == size(genes_start_end_TPM,1)
            intergenic = [intergenic; chromosome(end,1)-genes_start_end_TPM(end,2)];
        else
            intergenic = [intergenic; genes_start_end_TPM(j,1)-genes_start_end_TPM(j-1,2)];
        end
    end
    
    mean_intergenic_length_kawagutii = [mean_intergenic_length_kawagutii; mean(intergenic)];
    
    strand_changes = find(diff(sign(genes_start_end_TPM(:,3))))+1;

    tandem_gene_array_lengths = [];
    tandem_gene_array_lengths = [tandem_gene_array_lengths; strand_changes(1)];
    tandem_gene_array_lengths = [tandem_gene_array_lengths; diff(strand_changes)];

    mean_tandem_gene_array_length_kawagutii = [mean_tandem_gene_array_length_kawagutii; mean(tandem_gene_array_lengths)];

    strand_changes = [1; strand_changes];

    tandem_gene_array_lengths_bp = [];
    for j = 1:1:size(strand_changes,1)-1
        tandem_gene_array_lengths_bp = [tandem_gene_array_lengths_bp; genes_start_end_TPM(strand_changes(j+1)-1,2)-genes_start_end_TPM(strand_changes(j),1)];
    end

    mean_tandem_gene_array_length_bp_kawagutii = [mean_tandem_gene_array_length_bp_kawagutii; mean(tandem_gene_array_lengths_bp)];

    [coeff,score] = pca(chromosome(:,2:4));
    chromosome_PCA = score;

    %calculate moment of inertia tensor
    %Following: Chu, X., & Wang, J. (2023). Quantifying the large-scale chromosome structural dynamics during the mitosis-to-G1 phase transition of cell cycle. Open Biology, 13(11), 230175.
    Itensor=[[chromosome_PCA(:,1)'*chromosome_PCA(:,1) chromosome_PCA(:,1)'*chromosome_PCA(:,2) chromosome_PCA(:,1)'*chromosome_PCA(:,3)]; [chromosome_PCA(:,2)'*chromosome_PCA(:,1) chromosome_PCA(:,2)'*chromosome_PCA(:,2) chromosome_PCA(:,2)'*chromosome_PCA(:,3)]; [chromosome_PCA(:,3)'*chromosome_PCA(:,1) chromosome_PCA(:,3)'*chromosome_PCA(:,2) chromosome_PCA(:,3)'*chromosome_PCA(:,3)]];

    eigenvalues = eig(Itensor);

    asphericity_tensor_kawagutii = [asphericity_tensor_kawagutii; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];

    P_axis_1=[P_axis_1; prctile(abs(score(:,1)),95)];
    P_axis_2=[P_axis_2; prctile(abs(score(:,2)),95)];
    P_axis_3=[P_axis_3; prctile(abs(score(:,3)),95)];

    %get xyz coordinates for every base pair by linearly interpolating between the model
    chromosome_PCA_interpolated = zeros(chromosome(end,1),4);
    chromosome_PCA_interpolated(:,1) = linspace(1,chromosome(end,1),chromosome(end,1))';
    for j = 1:1:size(chromosome_PCA,1)-1
        chromosome_PCA_interpolated(HiC_resolution*(j-1)+1:HiC_resolution*j,2:4) = [linspace(chromosome_PCA(j,1),chromosome_PCA(j+1,1),HiC_resolution)' linspace(chromosome_PCA(j,2),chromosome_PCA(j+1,2),HiC_resolution)' linspace(chromosome_PCA(j,3),chromosome_PCA(j+1,3),HiC_resolution)'];
    end

    while genes_start_end_TPM(end,2) > chromosome_PCA_interpolated(end,1)
        genes_start_end_TPM(end,:)=[];
    end

    %genes_xyz =[];
    genes_center_xyz = [];

    %store coordinates for each base pair in each gene
    for k = 1:1:size(genes_start_end_TPM,1)
        %genes_xyz = [genes_xyz; chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),1) chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),2) chromosome_PCA_interpolated(genes_start_end_TPM(i,1):genes_start_end_TPM(i,2),3)];
        genes_center_xyz = [genes_center_xyz; mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),2)) mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),3)) mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),4))];
    end

    %column 1 is PC1, the chromosome long axis
    %cart2pol does x,y,z -> theta,rho,z. z is unchanged
    [theta_TPM, rho_TPM, z_TPM] = cart2pol(genes_center_xyz(:,2),genes_center_xyz(:,3),genes_center_xyz(:,1));

    %column 1 is PC1, the chromosome long axis
    [theta_chromosome, rho_chromosome, z_chromosome] = cart2pol(chromosome_PCA(:,2),chromosome_PCA(:,3),chromosome_PCA(:,1));

    temp_telomeres = [chromosome_PCA(1,:); chromosome_PCA(end,:)];

    %check that telomeres span a significant amount of clylindrical axis
    %total span of chromosome is [-1 1] (interval of length 2)
    %if  abs(z_chromosome(1,:)./max(abs(z_chromosome)) - z_chromosome(end,:)./max(abs(z_chromosome))) > 1
    if  abs(z_chromosome(1,:)./max(abs(z_chromosome)) - z_chromosome(end,:)./max(abs(z_chromosome))) > 0

        %rescale chromosome dimensions
        all_expression_agg = [all_expression_agg; rho_TPM./max(rho_chromosome) z_TPM./max(abs(z_chromosome)) abs(genes_start_end_TPM(:,3)) ones(size(genes_start_end_TPM,1),1).*i];
        telomeres = [telomeres; rho_chromosome(1,:)./max(rho_chromosome) z_chromosome(1,:)./max(abs(z_chromosome))];
        telomeres = [telomeres; rho_chromosome(end,:)./max(rho_chromosome) z_chromosome(end,:)./max(abs(z_chromosome))];
    end
    i
end

P=prctile(all_expression_agg(:,3),[20 80]);
%low expression <20th percentile
low_expression_agg = all_expression_agg(find(all_expression_agg(:,3)<P(1)),:);

%medium expression
med_expression_agg = all_expression_agg(find(all_expression_agg(:,3)>P(1) & all_expression_agg(:,3)<P(2)),:);

%high expression >80th percentile
high_expression_agg = all_expression_agg(find(all_expression_agg(:,3)>P(2)),:);

figure
[N,c] = hist3(low_expression_agg(:,1:2),[bins bins]);
N=flipud(N);
c{1}=fliplr(c{1});
for i = 1:1:size(N,1)
    N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
end
imagesc(N)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
hold on
%plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k')
hold off
%colormap("sky")
colormap(custom_cmap('green'))
cb = colorbar;
cb.Label.String = 'Normalized Gene Density';
cb.FontSize = 16;
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('cylidrical axis','FontSize', 24)
ylabel('radial axis','FontSize', 24)
t=get(cb,'Limits');
set(cb,'Ticks',[t(1),t(2)])
set(cb,'XTickLabel',{'Low','High',});
%xlim([2.5 20.5])

% figure
% [N,c] = hist3(med_expression_agg(:,1:2),[bins bins]);
% N=flipud(N);
% c{1}=fliplr(c{1});
% for i = 1:1:size(N,1)
%     N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
% end
% imagesc(N)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
% hold on
% plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k')
% hold off
% colormap("sky")
% cb = colorbar;
% cb.Label.String = 'Normalized Gene Density';
% cb.FontSize = 16;
% set(gca,'XTick',[])
% set(gca,'YTick',[])
% xlabel('cylidrical axis','FontSize', 24)
% ylabel('radial axis','FontSize', 24)
% t=get(cb,'Limits');
% set(cb,'Ticks',[t(1),t(2)])
% set(cb,'XTickLabel',{'Low','High',});

figure
[N,c] = hist3(high_expression_agg(:,1:2),[bins bins]);
N=flipud(N);
c{1}=fliplr(c{1});
for i = 1:1:size(N,1)
    N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
end
imagesc(N)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
hold on
%plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k')
hold off
%colormap("sky")
colormap(custom_cmap('green'))
cb = colorbar;
cb.Label.String = 'Normalized Gene Density';
cb.FontSize = 16;
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('cylidrical axis','FontSize', 24)
ylabel('radial axis','FontSize', 24)
t=get(cb,'Limits');
set(cb,'Ticks',[t(1),t(2)])
set(cb,'XTickLabel',{'Low','High',});
%xlim([2.5 20.5])

figure
[N,c] = hist3(all_expression_agg(:,1:2),[bins bins]);
N=flipud(N);
c{1}=fliplr(c{1});
for i = 1:1:size(N,1)
    N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
end
imagesc(N)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
hold on
plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k')
hold off
%colormap("sky")
colormap(custom_cmap('green'))
cb = colorbar;
cb.Label.String = 'Normalized Gene Density';
cb.FontSize = 16;
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('cylidrical axis','FontSize', 24)
ylabel('radial axis','FontSize', 24)
t=get(cb,'Limits');
set(cb,'Ticks',[t(1),t(2)])
set(cb,'XTickLabel',{'Low','High',});
%xlim([2.5 20.5])

ratio_P_axis_1_2_kawagutii = P_axis_1./P_axis_2;
ratio_P_axis_1_3_kawagutii = P_axis_1./P_axis_3;
ratio_P_axis_2_3_kawagutii = P_axis_2./P_axis_3;

for i = 1:1:num_chroms
    asphericity_kawagutii = [asphericity_kawagutii; max([ratio_P_axis_1_2_kawagutii(i), ratio_P_axis_1_3_kawagutii(i), ratio_P_axis_2_3_kawagutii(i)])];
end

all_expression_agg_kawagutii = all_expression_agg;

P_kawagutii=prctile(all_expression_agg_kawagutii(:,3),1:1:100);
P_microadriaticum=prctile(all_expression_agg_microadriaticum(:,3),1:1:100);

temp = [];
temp2 = [];
percentile_TPM_kawagutii = [];
for i = 1:1:num_chroms
    temp = all_expression_agg_kawagutii(find(all_expression_agg_kawagutii(:,4)==i),3);
    for j = 1:1:size(temp,1)
        temp2 = [temp2; mean(P_kawagutii<=temp(j))*100];
    end
    percentile_TPM_kawagutii = [percentile_TPM_kawagutii; mean(temp2)];
    temp2 = [];
end

temp = [];
temp2 = [];
percentile_TPM_microadriaticum = [];
for i = 1:1:num_chroms
    temp = all_expression_agg_microadriaticum(find(all_expression_agg_microadriaticum(:,4)==i),3);
    for j = 1:1:size(temp,1)
        temp2 = [temp2; mean(P_microadriaticum<=temp(j))*100];
    end
    percentile_TPM_microadriaticum = [percentile_TPM_microadriaticum; mean(temp2)];
    temp2 = [];
end

figure
hold on
%scatter(chromosome_length_microadriaticum./10^6, asphericity_microadriaticum, 25, 'blue', 'filled')
%scatter(chromosome_length_kawagutii./10^6, asphericity_kawagutii, 25, 'green', 'filled')
scatter(chromosome_length_microadriaticum./10^6, asphericity_tensor_microadriaticum, 25, 'blue', 'filled')
scatter(chromosome_length_kawagutii./10^6, asphericity_tensor_kawagutii, 25, 'green', 'filled')
lgd=legend({sprintf('Symbiodinium microadriaticum'),sprintf('Symbiodinium kawagutii')});
legend boxoff
hold off
ax = gca;
ax.FontSize = 16;
%ylim([1 2])
xlabel('Chromosome Length [Mbp]','FontSize', 24)
ylabel('Chromosome Asphericity','FontSize', 24)

figure
hold on
%scatter(median_TPM_microadriaticum, asphericity_microadriaticum, 25, 'blue', 'filled')
%scatter(median_TPM_kawagutii, asphericity_kawagutii, 25, 'green', 'filled')
scatter(median_TPM_microadriaticum, asphericity_tensor_microadriaticum, 25, 'blue', 'filled')
scatter(median_TPM_kawagutii, asphericity_tensor_kawagutii, 25, 'green', 'filled')
lgd=legend({sprintf('Symbiodinium microadriaticum'),sprintf('Symbiodinium kawagutii')});
legend boxoff
hold off
ax = gca;
ax.FontSize = 16;
%ylim([1 2])
xlabel('Mean TPM','FontSize', 24)
ylabel('Chromosome Asphericity','FontSize', 24)

figure
hold on
%scatter(percentile_TPM_microadriaticum, asphericity_microadriaticum, 25, 'blue', 'filled')
%scatter(percentile_TPM_kawagutii, asphericity_kawagutii, 25, 'green', 'filled')
scatter(percentile_TPM_microadriaticum, asphericity_tensor_microadriaticum, 25, 'blue', 'filled')
scatter(percentile_TPM_kawagutii, asphericity_tensor_kawagutii, 25, 'green', 'filled')
lgd=legend({sprintf('Symbiodinium microadriaticum'),sprintf('Symbiodinium kawagutii')});
legend boxoff
hold off
ax = gca;
ax.FontSize = 16;
%ylim([1 2])
xlabel('Average TPM Percentile Per Chromosome','FontSize', 24)
ylabel('Chromosome Asphericity','FontSize', 24)

figure
hold on
%scatter(percent_gene_content_microadriaticum, asphericity_microadriaticum, 25, 'blue', 'filled')
%scatter(percent_gene_content_kawagutii, asphericity_kawagutii, 25, 'green', 'filled')
scatter(percent_gene_content_microadriaticum, asphericity_tensor_microadriaticum, 25, 'blue', 'filled')
scatter(percent_gene_content_kawagutii, asphericity_tensor_kawagutii, 25, 'green', 'filled')
lgd=legend({sprintf('Symbiodinium microadriaticum'),sprintf('Symbiodinium kawagutii')});
legend boxoff
hold off
ax = gca;
ax.FontSize = 16;
%ylim([1 2])
xlabel('% Gene Content','FontSize', 24)
ylabel('Chromosome Asphericity','FontSize', 24)

figure
hold on
%scatter(mean_intergenic_length_microadriaticum./10^3, asphericity_microadriaticum, 25, 'blue', 'filled')
%scatter(mean_intergenic_length_kawagutii./10^3, asphericity_kawagutii, 25, 'green', 'filled')
scatter(mean_intergenic_length_microadriaticum./10^3, asphericity_tensor_microadriaticum, 25, 'blue', 'filled')
scatter(mean_intergenic_length_kawagutii./10^3, asphericity_tensor_kawagutii, 25, 'green', 'filled')
lgd=legend({sprintf('Symbiodinium microadriaticum'),sprintf('Symbiodinium kawagutii')});
legend boxoff
hold off
ax = gca;
ax.FontSize = 16;
%ylim([1 2])
xlabel('Mean intergenic length [kbp]','FontSize', 24)
ylabel('Chromosome Asphericity','FontSize', 24)

figure
hold on
%scatter(mean_tandem_gene_array_length_microadriaticum, asphericity_microadriaticum, 25, 'blue', 'filled')
%scatter(mean_tandem_gene_array_length_kawagutii, asphericity_kawagutii, 25, 'green', 'filled')
scatter(mean_tandem_gene_array_length_microadriaticum, asphericity_tensor_microadriaticum, 25, 'blue', 'filled')
scatter(mean_tandem_gene_array_length_kawagutii, asphericity_tensor_kawagutii, 25, 'green', 'filled')
lgd=legend({sprintf('Symbiodinium microadriaticum'),sprintf('Symbiodinium kawagutii')});
legend boxoff
hold off
ax = gca;
ax.FontSize = 16;
%ylim([1 2])
xlabel('Mean Tandem Gene Array Length','FontSize', 24)
ylabel('Chromosome Asphericity','FontSize', 24)

figure
hold on
%scatter(mean_tandem_gene_array_length_bp_microadriaticum./10^3, asphericity_microadriaticum, 25, 'blue', 'filled')
%scatter(mean_tandem_gene_array_length_bp_kawagutii./10^3, asphericity_kawagutii, 25, 'green', 'filled')
scatter(mean_tandem_gene_array_length_bp_microadriaticum./10^3, asphericity_tensor_microadriaticum, 25, 'blue', 'filled')
scatter(mean_tandem_gene_array_length_bp_kawagutii./10^3, asphericity_tensor_kawagutii, 25, 'green', 'filled')
lgd=legend({sprintf('Symbiodinium microadriaticum'),sprintf('Symbiodinium kawagutii')});
legend boxoff
hold off
ax = gca;
ax.FontSize = 16;
%ylim([1 2])
xlabel('Mean Tandem Gene Array Length [kbp]','FontSize', 24)
ylabel('Chromosome Asphericity','FontSize', 24)

mean(all_expression_agg(find(all_expression_agg(:,3)>P(2)),3))/mean(all_expression_agg(find(all_expression_agg(:,3)<P(1)),3)) %(mean (high TPM value))/(mean (low TPM value))

figure
hold on
[counts1, binCenters1] = hist(ratio_P_axis_1_2_microadriaticum, 10);
[counts2, binCenters2] = hist(ratio_P_axis_1_3_microadriaticum, 10);
[counts3, binCenters3] = hist(ratio_P_axis_2_3_microadriaticum, 10);
[counts4, binCenters4] = hist(ratio_P_axis_1_2_kawagutii, 10);
[counts5, binCenters5] = hist(ratio_P_axis_1_3_kawagutii, 10);
[counts6, binCenters6] = hist(ratio_P_axis_2_3_kawagutii, 10);
plot(binCenters1, counts1, 'g-', 'LineWidth',2);
hold on;
plot(binCenters2, counts2, 'g-','LineWidth',2);
plot(binCenters3, counts3, 'g-','LineWidth',2);
plot(binCenters4, counts4, 'b-','LineWidth',2);
plot(binCenters5, counts5, 'b-','LineWidth',2);
plot(binCenters6, counts6, 'b-','LineWidth',2);
xlabel('Asphericity (Ratio of Principle Components)','FontSize', 24)
ylabel('Count','FontSize', 24)
ax = gca;
ax.FontSize = 16;
xlim([1 2])
% Put up legend.
legend1 = sprintf('Symbiodinium kawagutii');
legend4 = sprintf('Symbiodinium microadriaticum');
lgd=legend({legend1, '', '', legend4, '', ''});
legend boxoff
lgd.FontSize = 20;

figure
hold on
[counts1, binCenters1] = hist(asphericity_tensor_kawagutii, 10);
[counts2, binCenters2] = hist(asphericity_tensor_microadriaticum, 10);
plot(binCenters1, counts1, 'g-', 'LineWidth',2);
plot(binCenters2, counts2, 'b-','LineWidth',2);
xlabel('Asphericity','FontSize', 24)
ylabel('Count','FontSize', 24)
ax = gca;
ax.FontSize = 16;
xlim([0 0.25])
% Put up legend.
legend1 = sprintf('Symbiodinium kawagutii');
legend2 = sprintf('Symbiodinium microadriaticum');
lgd=legend({legend1, legend2});
legend boxoff
lgd.FontSize = 20;

% %genes per unit volume plot
% chrom_xyz = chromosome_PCA_interpolated;
%
% %normalized gene density
% %% Compute density
% % Put points into 3D bins; xyzBinNum is an nx3 matrix containing
% % the bin ID for n values in xyz for the [x,y,z] axes.
% nBins = 10;  % number of bins
% chrom_xbins = linspace(min(chrom_xyz(:,1)),max(chrom_xyz(:,1))*1,nBins+1);
% chrom_ybins = linspace(min(chrom_xyz(:,2)),max(chrom_xyz(:,2))*1,nBins+1);
% chrom_zbins = linspace(min(chrom_xyz(:,3)),max(chrom_xyz(:,3))*1,nBins+1);
% chrom_xyzBinNum = [...
%     discretize(chrom_xyz(:,1),chrom_xbins), ...
%     discretize(chrom_xyz(:,2),chrom_ybins), ...
%     discretize(chrom_xyz(:,3),chrom_zbins), ...
%     ];
%
% genes_xyzBinNum = [...
%     discretize(genes_xyz(:,1),chrom_xbins), ...
%     discretize(genes_xyz(:,2),chrom_ybins), ...
%     discretize(genes_xyz(:,3),chrom_zbins), ...
%     ];
% % bin3D is a mx3 matrix of m unique 3D bins that appear
% % in xyzBinNum, sorted.  binNum is a nx1 vector of bin
% % numbers identifying the bin for each xyz point. For example,
% % b=xyz(j,:) belongs to bins3D(b,:).
% [chrom_bins3D, ~, chrom_binNum] = unique(chrom_xyzBinNum, 'rows');
%
% %map gene_density onto the same bins as chrom_density
% gene_density = [];
% for i = 1:1:size(chrom_bins3D,1)
%         gene_density = [gene_density; sum(ismember(genes_xyzBinNum, chrom_bins3D(i,:),'rows'));];
%         size(chrom_bins3D,1)-i
% end
%
% % density is a mx1 vector of integers showing the number of
% % xyz points in each of the bins3D. To see the number of points
% % in bins3D(k,:), density(k).
% chrom_density = histcounts(chrom_binNum,[1:size(chrom_bins3D,1),inf])';
%
% % Compute bin centers
% chrom_xbinCnt = chrom_xbins(2:end)-diff(chrom_xbins)/2;
% chrom_ybinCnt = chrom_ybins(2:end)-diff(chrom_ybins)/2;
% chrom_zbinCnt = chrom_zbins(2:end)-diff(chrom_zbins)/2;
%
% chrom_bins3D(find(gene_density==0),:) = [];
% chrom_density(find(gene_density==0)) = [];
% gene_density(find(gene_density==0)) = [];
%
% %% Plot raw data
% fig = figure();
% %% Plot scatter3
% scatter3(...
%     chrom_xbinCnt(chrom_bins3D(:,1)), ...
%     chrom_ybinCnt(chrom_bins3D(:,2)), ...
%     chrom_zbinCnt(chrom_bins3D(:,3)), ...
%     gene_density./chrom_density*100, ...
%     gene_density./chrom_density, 'filled', ...
%     'MarkerFaceAlpha',.6)
% axis equal
% box on
% xlabel('x')
% ylabel('y')
% zlabel('z')
% cb = colorbar
% cb.Label.String = 'Gene Density/Chromsome Density';
% cb.FontSize = 16;


%see: https://github.com/JonathanRob/GeneSetAnalysisMatlab/blob/master/custom_cmap.m
function cmap = custom_cmap(colors,n,fractions)
%custom_cmap  Generate a custom or pre-defined RGB colormap.
%
% Usage:
%
%   cmap = custom_cmap(colors,n,fractions);
%
% Input:
%
%   colors      A Cx3 matrix of C RGB values to include in the colormap, or
%               a string specifying a pre-built custom colormap:
%               'redblue', 'redbluebright', 'green', 'blue', 'darkblue',
%               'magma'
%
%               ** custom_cmap('preview'); shows a preview of the pre-built
%                  colormap options
%
%   n           Total number of colors desired in the final colormap.
%               Default = 100.
%
%   fractions   Numeric vector representing the fractions of the colormap
%               dedicated to each color (should sum to 1). If empty or not
%               provided, the colors will be distributed evenly.
%
% Output:
%
%   cmap        An Nx3 matrix of RGB values comprising a colormap.
%
%
% Jonathan Robinson, 2020-02-03

if nargin < 2 || isempty(n)
    n = 100;
end

if ~isnumeric(colors)

    % define custom colormaps
    switch lower(colors)
        case 'redblue'
            hotcmap = hot(round(60*n/100));
            cmap = flipud([hotcmap(round(10*n/100+1):end,:); custom_cmap([0.1 0.1 0.8;0 0.7 0.9;1 1 1],round(n/2),[0.25;0.4;0.35])]);
        case 'redbluebright'
            hotcmap = hot(round(80*n/100));
            cmap = flipud([hotcmap(round(30*n/100+1):end,:); custom_cmap([0.2 0.2 1;0 1 1;1 1 1],round(n/2),[0.25;0.4;0.35])]);
        case 'green'
            cmap = custom_cmap([0.0118 0.2392 0.0157;0.2392 0.7490 0.2471;1 1 1],n,[0.25;0.4;0.35]);
        case 'blue'
            cmap = custom_cmap([8,58,106;8,81,156;49,130,189;107,174,214;189,215,231;255,255,255]./255,n);
        case 'darkblue'
            cmap = custom_cmap([0.015,0.106,0.247; 0.031,0.227,0.416; 0.055,0.361,0.580; 0.400,0.592,0.690; 0.706,0.804,0.875],n);
        case 'magma'
            cmap = custom_cmap([15,15,19; 42,0,92; 117,15,110; 211,51,86; 253,144,91; 251,255,178]./255,n);
        case 'whitemagma'
            cmap = custom_cmap([15,15,19; 42,0,92; 117,15,110; 211,51,86; 253,144,91; 251,255,178; 255,255,255]./255,n,...
                [0.165;0.165;0.165;0.165;0.165;0.165;0.01]);
            cmap(1,:) = [1,1,1];  % ensure that last value is pure white
        case 'redbluedark'
            cmap = custom_cmap([0.573,0.216,0.212; 1,1,1; 0.039,0.345,0.529],n);
        case 'preview'
            % preview all pre-defined colormap options
            cmap_list = {'redblue';'redbluebright';'redbluedark';'green';'blue';'darkblue';'magma';'whitemagma'};
            Ncmap = length(cmap_list);
            cmap = [];
            for i = 1:length(cmap_list)
                cmap = [cmap; custom_cmap(cmap_list{i},200)];
            end
            plotmat = flipud(reshape(0:(200*Ncmap-1),200,length(cmap_list))');
            plotmat(end+1,end+1) = 0;
            figure;
            h = pcolor(plotmat); hold on
            set(h,'EdgeColor','none');
            set(gca,'YTick',(1:Ncmap)+0.5,'YTickLabels',flipud(cmap_list),'FontSize',14);
            plot([zeros(1,Ncmap-1);202*ones(1,Ncmap-1)],repmat(2:Ncmap,2,1),'k-');
            colormap(cmap);
            return
        otherwise
            error('"%s" is not a valid colormap option.',colors);
    end

else

    [numcolors,z] = size(colors);
    if z ~= 3
        error('COLORS must be an N x 3 matrix.');
    end

    if nargin < 3 || isempty(fractions) || length(unique(fractions)) == 1
        cmap = flipud(interp1(1:numcolors, colors, linspace(1,numcolors,n)));
    else
        if n < numcolors
            warning('Too many colors. Not all will be included in CMAP.');
        elseif length(fractions) ~= numcolors
            error('COLORS and FRACTIONS must contain equal number of elements');
        end

        cmap = zeros(n,3);
        cx = linspace(0,1,n+1)' + 1/(2*n);
        cx(end) = [];

        x = sort([0; cumsum(fractions); cumsum(fractions(1:end-2)) + fractions(2:end-1)/2 ]);
        y = zeros(numcolors*2 - 1,3);
        for i = 1:numcolors
            y(i*2-1,:) = colors(i,:);
            if i < numcolors
                y(i*2,:) = mean([colors(i,:);colors(i+1,:)]);
            end
        end
        for i = 1:3
            cmap(:,i) = interp1(x,y(:,i),cx);
        end
        cmap = flipud(cmap);
    end
end
end


