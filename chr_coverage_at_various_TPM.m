%Fig S8
%Cumulative chromosomal coverage of active gene sequence with TPM > x. 
%chromosome surface area to volume ratio

num_chroms = 50;
x = logspace(-3,4);
expression_coverage = zeros(length(x),num_chroms);
exposed_shell_deep = [];
exposed_shell_med = [];
exposed_shell_shallow = [];

%microadriaticum
for i=1:num_chroms

    chromosome = importdata(sprintf('symbiodinium_microadriaticum_chr%i_3D.xyz',i));
    RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));

    chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
    chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));
    chromosome(:,4) = chromosome(:,4) - mean(chromosome(:,4));

    TPM_temp = RNAseq.data(:,8);
    RNAseq_temp = RNAseq.data;
    for j = 1:length(x)
        expression_coverage(j,i) = sum(RNAseq_temp(find(abs(TPM_temp)>x(j)),2)-RNAseq_temp(find(abs(TPM_temp)>x(j)),1))/(size(chromosome,1)*5000);
    end

    [coeff,score] = pca(chromosome(:,2:4));
    chromosome_PCA = score;

    P_axis_1=prctile(chromosome_PCA(:,1)-min(chromosome_PCA(:,1)),95);
    P_axis_2=prctile(chromosome_PCA(:,2)-min(chromosome_PCA(:,2)),95);
    P_axis_3=prctile(chromosome_PCA(:,3)-min(chromosome_PCA(:,3)),95);

    volume = 4*pi*(P_axis_1/2)*(P_axis_2/2)*(P_axis_3/2)/3;
    volume_core = 4*pi*((P_axis_1/2*(1-0.3))*(P_axis_2/2*(1-0.3))*(P_axis_3/2*(1-0.3)))/3;
    exposed_shell_deep = [exposed_shell_deep; ((volume-volume_core)/volume)*expression_coverage(1,i)]; %not all DNA is a gene, multiply by active gene density

    volume_core = 4*pi*((P_axis_1/2*(1-0.2))*(P_axis_2/2*(1-0.2))*(P_axis_3/2*(1-0.2)))/3;
    exposed_shell_med = [exposed_shell_med; ((volume-volume_core)/volume)*expression_coverage(1,i)]; %not all DNA is a gene, multiply by active gene density

    volume_core = 4*pi*((P_axis_1/2*(1-0.1))*(P_axis_2/2*(1-0.1))*(P_axis_3/2*(1-0.1)))/3;
    exposed_shell_shallow = [exposed_shell_shallow; ((volume-volume_core)/volume)*expression_coverage(1,i)]; %not all DNA is a gene, multiply by active gene density
end

figure
for i = 1:num_chroms
    semilogx(x,expression_coverage(:,i).*100, 'Color', [0 0.4470 0.7410]) %multiplication by 100 to convert fraction to percentage
    yline(exposed_shell_deep(i).*100,'--r')
    yline(exposed_shell_med(i).*100,'--y')
    yline(exposed_shell_shallow(i).*100,'--g')
    hold on
end
hxl = yline(exposed_shell_deep(end).*100,'--k','30% Radial Penetration')
hxl.FontSize = 14;
hxl = yline(exposed_shell_med(end).*100,'--k','20% Radial Penetration')
hxl.FontSize = 14;
hxl = yline(exposed_shell_shallow(end).*100,'--k','10% Radial Penetration')
hxl.FontSize = 14;
ylabel('% Chromosome Sequence with TPM > x')
xlabel('x [TPM]')
ax = gca;
ax.FontSize = 16;
xlim([10^(-3) 10^4])
ylim([0 100])

x = logspace(-3,4);
expression_coverage = zeros(length(x),num_chroms);
exposed_shell_deep = [];
exposed_shell_med = [];
exposed_shell_shallow = [];

for i=1:num_chroms

    chromosome = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i.xyz',i));
    RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));

    chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
    chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));
    chromosome(:,4) = chromosome(:,4) - mean(chromosome(:,4));

    TPM_temp = RNAseq.data(:,8);
    RNAseq_temp = RNAseq.data;
    for j = 1:length(x)
        expression_coverage(j,i) = sum(RNAseq_temp(find(abs(TPM_temp)>x(j)),2)-RNAseq_temp(find(abs(TPM_temp)>x(j)),1))/(size(chromosome,1)*5000);
    end

    [coeff,score] = pca(chromosome(:,2:4));
    chromosome_PCA = score;

    P_axis_1=prctile(chromosome_PCA(:,1)-min(chromosome_PCA(:,1)),95);
    P_axis_2=prctile(chromosome_PCA(:,2)-min(chromosome_PCA(:,2)),95);
    P_axis_3=prctile(chromosome_PCA(:,3)-min(chromosome_PCA(:,3)),95);

    volume = 4*pi*(P_axis_1/2)*(P_axis_2/2)*(P_axis_3/2)/3;
    volume_core = 4*pi*((P_axis_1/2*(1-0.3))*(P_axis_2/2*(1-0.3))*(P_axis_3/2*(1-0.3)))/3;
    exposed_shell_deep = [exposed_shell_deep; ((volume-volume_core)/volume)*expression_coverage(1,i)]; %not all DNA is a gene, multiply by active gene density

    volume_core = 4*pi*((P_axis_1/2*(1-0.2))*(P_axis_2/2*(1-0.2))*(P_axis_3/2*(1-0.2)))/3;
    exposed_shell_med = [exposed_shell_med; ((volume-volume_core)/volume)*expression_coverage(1,i)]; %not all DNA is a gene, multiply by active gene density

    volume_core = 4*pi*((P_axis_1/2*(1-0.1))*(P_axis_2/2*(1-0.1))*(P_axis_3/2*(1-0.1)))/3;
    exposed_shell_shallow = [exposed_shell_shallow; ((volume-volume_core)/volume)*expression_coverage(1,i)]; %not all DNA is a gene, multiply by active gene density
end

%kawagutii
figure
for i = 1:num_chroms
    semilogx(x,expression_coverage(:,i).*100, 'Color', [0.4660 0.6740 0.1880]) %multiplication by 100 to convert fraction to percentage
    yline(exposed_shell_deep(i).*100,'--r')
    yline(exposed_shell_med(i).*100,'--y')
    yline(exposed_shell_shallow(i).*100,'--g')
    hold on
end
hxl = yline(exposed_shell_deep(end).*100,'--k','30% Radial Penetration')
hxl.FontSize = 14;
hxl = yline(exposed_shell_med(end).*100,'--k','20% Radial Penetration')
hxl.FontSize = 14;
hxl = yline(exposed_shell_shallow(end).*100,'--k','10% Radial Penetration')
hxl.FontSize = 14;
ylabel('% Chromosome Sequence with TPM > x')
xlabel('x [TPM]')
ax = gca;
ax.FontSize = 16;
xlim([10^(-3) 10^4])
ylim([0 100])