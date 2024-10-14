%Fig S10
%Cumulative chromosomal coverage of active gene sequence with TPM > x. 
%chromosome surface area to volume ratio

num_chroms = 75;
x = logspace(-3,4);
expression_coverage = zeros(length(x),num_chroms);
SA_volume_ratio = [];

%microadriaticum
for i=1:num_chroms

    chromosome = importdata(sprintf('symbiodinium_microadriaticum_chr%i_3D.xyz',i));
    RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));

    chromosome(:,1) = chromosome(:,1) - mean(chromosome(:,1));
    chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
    chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));

    TPM_temp = RNAseq.data(:,8);
    RNAseq_temp = RNAseq.data;
    for j = 1:length(x)
        expression_coverage(j,i) = sum(RNAseq_temp(find(abs(TPM_temp)>x(j)),2)-RNAseq_temp(find(abs(TPM_temp)>x(j)),1))/(size(chromosome,1)*5000);
    end

    [coeff,score] = pca(chromosome);
    chromosome_PCA = score;

    P_axis_1=prctile(chromosome_PCA(:,1)-min(chromosome_PCA(:,1)),95);
    P_axis_2=prctile(chromosome_PCA(:,2)-min(chromosome_PCA(:,2)),95);
    P_axis_3=prctile(chromosome_PCA(:,3)-min(chromosome_PCA(:,3)),95);

    SA = 4*pi*(((P_axis_1*P_axis_2)^1.6+(P_axis_1*P_axis_3)^1.6+(P_axis_2*P_axis_3)^1.6)/3)^(1/1.6);
    volume = 4*pi*P_axis_1*P_axis_2*P_axis_3/3;
    SA_volume_ratio = [SA_volume_ratio; SA/volume];
end

figure
for i = 1:num_chroms
    semilogx(x,expression_coverage(:,i), 'Color', [0 0.4470 0.7410])
    yline(SA_volume_ratio(i),'--k')
    hold on
end
hxl = yline(SA_volume_ratio(end),'--k','Chr Surface Area/Volume')
hxl.FontSize = 16;
hxl.LabelHorizontalAlignment='left';
ylabel('% Chromosome Sequence with TPM > x')
xlabel('x [TPM]')
ax = gca;
ax.FontSize = 16;
xlim([10^(-3) 10^4])
ylim([0 1])

x = logspace(-3,4);
expression_coverage = zeros(length(x),num_chroms);
SA_volume_ratio = [];

for i=1:num_chroms

    chromosome = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i.xyz',i));
    RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));

    chromosome(:,1) = chromosome(:,1) - mean(chromosome(:,1));
    chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
    chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));

    TPM_temp = RNAseq.data(:,8);
    RNAseq_temp = RNAseq.data;
    for j = 1:length(x)
        expression_coverage(j,i) = sum(RNAseq_temp(find(abs(TPM_temp)>x(j)),2)-RNAseq_temp(find(abs(TPM_temp)>x(j)),1))/(size(chromosome,1)*5000);
    end

    [coeff,score] = pca(chromosome);
    chromosome_PCA = score;

    P_axis_1=prctile(chromosome_PCA(:,1)-min(chromosome_PCA(:,1)),95);
    P_axis_2=prctile(chromosome_PCA(:,2)-min(chromosome_PCA(:,2)),95);
    P_axis_3=prctile(chromosome_PCA(:,3)-min(chromosome_PCA(:,3)),95);

    SA = 4*pi*(((P_axis_1*P_axis_2)^1.6+(P_axis_1*P_axis_3)^1.6+(P_axis_2*P_axis_3)^1.6)/3)^(1/1.6);
    volume = 4*pi*P_axis_1*P_axis_2*P_axis_3/3;
    SA_volume_ratio = [SA_volume_ratio; SA/volume];
end

%kawagutii
figure
for i = 1:num_chroms
    semilogx(x,expression_coverage(:,i), 'Color', [0.4660 0.6740 0.1880])
    yline(SA_volume_ratio(i),'--k')
    hold on
end
hxl = yline(SA_volume_ratio(end),'--k','Chr Surface Area/Volume')
hxl.FontSize = 16;
hxl.LabelHorizontalAlignment='left';
ylabel('% Chromosome Sequence with TPM > x')
xlabel('x [TPM]')
ax = gca;
ax.FontSize = 16;
xlim([10^(-3) 10^4])
ylim([0 1])