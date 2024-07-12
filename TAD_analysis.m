clc
clear

num_chroms = 52;
dist_to_TAD = 45000;

asphericity_tensor_TADs_microadriaticum =[];
asphericity_tensor_microadriaticum = [];

asphericity_tensor_TADs_kawagutii=[];
asphericity_tensor_kawagutii = [];

TAD_TPM_microadriaticum = [];
TAD_TPM_kawagutii = [];

TAD_num_microadriaticum = [];
TAD_percent_coverage_microadriaticum = [];

TAD_num_kawagutii = [];
TAD_percent_coverage_kawagutii = [];

P_axis_1_microadriaticum = [];
P_axis_2_microadriaticum = [];
P_axis_3_microadriaticum = [];

P_axis_1_kawagutii = [];
P_axis_2_kawagutii = [];
P_axis_3_kawagutii = [];

for i=1:num_chroms
    TAD_boundary_microadriaticum = [];
    TADs_microadriaticum = [];
    temp_TAD_boundary = [];

    if isfile(sprintf('symbiodinium_microadriaticum_chr%i.txt.bed',i))
        % File exists.
        chromosome = importdata(sprintf('symbiodinium_microadriaticum_chr%i_3D.xyz',i));
        RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));
        TADfile = importdata(sprintf('symbiodinium_microadriaticum_chr%i.txt.bed',i));

        chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
        chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));
        chromosome(:,4) = chromosome(:,4) - mean(chromosome(:,4));

        HiC_resolution = chromosome(2,1)-chromosome(1,1);

        genes_start_end_TPM = [RNAseq.data(:,1) RNAseq.data(:,2) RNAseq.data(:,8)];

        for j = 1:1:size(TADfile,1)
            istab=strfind(TADfile{j},char(9));
            TAD_boundary_microadriaticum = [TAD_boundary_microadriaticum; str2double(extractBetween(TADfile{j},istab(1)+1,istab(2)-1)); str2double(extractBetween(TADfile{j},istab(2)+1,istab(3)-1))];
            temp_TAD_boundary=[temp_TAD_boundary; str2double(extractBetween(TADfile{j},istab(1)+1,istab(2)-1)) str2double(extractBetween(TADfile{j},istab(2)+1,istab(3)-1))];
        end

        TADs_microadriaticum = [TADs_microadriaticum; temp_TAD_boundary];

        TAD_num_microadriaticum = [TAD_num_microadriaticum; size(TADfile,1)];
        TAD_percent_coverage_microadriaticum = [TAD_percent_coverage_microadriaticum; (sum(temp_TAD_boundary(:,2)-temp_TAD_boundary(:,1))./chromosome(end,1)).*100];

        TPM_bp = zeros(chromosome(end,1),1);
        for j = 1:1:size(genes_start_end_TPM,1)
            TPM_bp(genes_start_end_TPM(j,1):genes_start_end_TPM(j,2)) = genes_start_end_TPM(j,3);
        end

        %100bp sliding windows
        mov_avg_TPM = movmedian(TPM_bp,100);

        for j = 1:1:size(TAD_boundary_microadriaticum,1)
            %30000bp upstream and downstream of TAD
            if TAD_boundary_microadriaticum(j)-dist_to_TAD > 0 & TAD_boundary_microadriaticum(j)+dist_to_TAD < chromosome(end,1)

                left_avg_TPM = median(mov_avg_TPM(TAD_boundary_microadriaticum(j)-dist_to_TAD:TAD_boundary_microadriaticum(j))');
                right_avg_TPM = median(mov_avg_TPM(TAD_boundary_microadriaticum(j):TAD_boundary_microadriaticum(j)+dist_to_TAD)');

                if right_avg_TPM > left_avg_TPM
                    TAD_TPM_microadriaticum = [TAD_TPM_microadriaticum; -1.*(mov_avg_TPM(TAD_boundary_microadriaticum(j)-dist_to_TAD:TAD_boundary_microadriaticum(j)+dist_to_TAD)')];
                else
                    TAD_TPM_microadriaticum = [TAD_TPM_microadriaticum; mov_avg_TPM(TAD_boundary_microadriaticum(j)-dist_to_TAD:TAD_boundary_microadriaticum(j)+dist_to_TAD)'];
                end
                %elseif TAD_boundary_microadriaticum(j)-dist_to_TAD < 0
                %TAD_TPM = [TAD_TPM; mov_avg_TPM(1:TAD_boundary_microadriaticum(j)+dist_to_TAD)'];
                %elseif TAD_boundary_microadriaticum(j)+dist_to_TAD > chromosome(end,1)
                %TAD_TPM = [TAD_TPM; mov_avg_TPM(TAD_boundary_microadriaticum(j)-dist_to_TAD:chromosome(end,1))'];
            end
        end

        [coeff,score] = pca(chromosome(:,2:4));
        chromosome_PCA = score;

        P_axis_1_microadriaticum=[P_axis_1_microadriaticum; prctile(abs(score(:,1)),95)];
        P_axis_2_microadriaticum=[P_axis_2_microadriaticum; prctile(abs(score(:,2)),95)];
        P_axis_3_microadriaticum=[P_axis_3_microadriaticum; prctile(abs(score(:,3)),95)];

        for j = 1:1:size(TADs_microadriaticum,1)
            TAD_PCA = chromosome_PCA(round(TADs_microadriaticum(j,1)./HiC_resolution)+1:round(TADs_microadriaticum(j,2)./HiC_resolution),:);
            %calculate moment of inertia tensor
            %Following: Chu, X., & Wang, J. (2023). Quantifying the large-scale chromosome structural dynamics during the mitosis-to-G1 phase transition of cell cycle. Open Biology, 13(11), 230175.
            Itensor_TAD = [[TAD_PCA(:,1)'*TAD_PCA(:,1) TAD_PCA(:,1)'*TAD_PCA(:,2) TAD_PCA(:,1)'*TAD_PCA(:,3)]; [TAD_PCA(:,2)'*TAD_PCA(:,1) TAD_PCA(:,2)'*TAD_PCA(:,2) TAD_PCA(:,2)'*TAD_PCA(:,3)]; [TAD_PCA(:,3)'*TAD_PCA(:,1) TAD_PCA(:,3)'*TAD_PCA(:,2) TAD_PCA(:,3)'*TAD_PCA(:,3)]];
            eigenvalues = eig(Itensor_TAD);
            asphericity_tensor_TADs_microadriaticum = [asphericity_tensor_TADs_microadriaticum; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];
        end
        Itensor_chromosome = [[chromosome_PCA(:,1)'*chromosome_PCA(:,1) chromosome_PCA(:,1)'*chromosome_PCA(:,2) chromosome_PCA(:,1)'*chromosome_PCA(:,3)]; [chromosome_PCA(:,2)'*chromosome_PCA(:,1) chromosome_PCA(:,2)'*chromosome_PCA(:,2) chromosome_PCA(:,2)'*chromosome_PCA(:,3)]; [chromosome_PCA(:,3)'*chromosome_PCA(:,1) chromosome_PCA(:,3)'*chromosome_PCA(:,2) chromosome_PCA(:,3)'*chromosome_PCA(:,3)]];
        eigenvalues = eig(Itensor_chromosome);
        asphericity_tensor_microadriaticum = [asphericity_tensor_microadriaticum; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];
        
        i
    else
        % File does not exist.
    end
end

%TPM versus distance from TAD boundary
figure
hold on
plot((-dist_to_TAD:1:dist_to_TAD)./1000,mean(TAD_TPM_microadriaticum,1), 'LineWidth', 8, 'Color', [0 0.4470 0.7410])
xline(0)
yline(0)
ylim([-8 6])
ax = gca;
ax.FontSize = 16;
xlabel('Distance from TAD Boundary [kbp]','FontSize', 24)
ylabel('TPM','FontSize', 24)

%TAD asphericity
figure
hold on
histogram(asphericity_tensor_TADs_microadriaticum, 20, FaceColor = [0 0.4470 0.7410]./3) %a lighter colour
histogram(asphericity_tensor_microadriaticum, FaceColor = [0 0.4470 0.7410])
lgd=legend({'TADs Only','Entire Chromosome'});
legend boxoff
hold off
xlim([0 1])
ax = gca;
ax.FontSize = 16;
xlabel('Asphericity','FontSize', 24)
ylabel('Count','FontSize', 24)

for i=1:num_chroms
    TAD_boundary_kawagutii = [];
    TADs_kawagutii = [];
    temp_TAD_boundary = [];

    if isfile(sprintf('s_kawagutii_V3_HiC_scaffold_%i_for_CSynth.txt.bed',i))
        % File exists.
        chromosome = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i.xyz',i));
        RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));
        TADfile = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i_for_CSynth.txt.bed',i));

        chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
        chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));
        chromosome(:,4) = chromosome(:,4) - mean(chromosome(:,4));

        HiC_resolution = chromosome(2,1)-chromosome(1,1);

        genes_start_end_TPM = [RNAseq.data(:,1) RNAseq.data(:,2) RNAseq.data(:,8)];

        for j = 1:1:size(TADfile,1)
            istab=strfind(TADfile{j},char(9));
            TAD_boundary_kawagutii = [TAD_boundary_kawagutii; str2double(extractBetween(TADfile{j},istab(1)+1,istab(2)-1)); str2double(extractBetween(TADfile{j},istab(2)+1,istab(3)-1))];
            temp_TAD_boundary=[temp_TAD_boundary; str2double(extractBetween(TADfile{j},istab(1)+1,istab(2)-1)) str2double(extractBetween(TADfile{j},istab(2)+1,istab(3)-1))];
        end

        TADs_kawagutii = [TADs_kawagutii; temp_TAD_boundary];

        TAD_num_kawagutii = [TAD_num_kawagutii; size(TADfile,1)];
        TAD_percent_coverage_kawagutii = [TAD_percent_coverage_kawagutii; (sum(temp_TAD_boundary(:,2)-temp_TAD_boundary(:,1))./chromosome(end,1)).*100];

        TPM_bp = zeros(chromosome(end,1),1);
        for j = 1:1:size(genes_start_end_TPM,1)
            TPM_bp(genes_start_end_TPM(j,1):genes_start_end_TPM(j,2)) = genes_start_end_TPM(j,3);
        end

        %100bp sliding windows
        mov_avg_TPM = movmedian(TPM_bp,100);

        for j = 1:1:size(TAD_boundary_kawagutii,1)
            %30000bp upstream and downstream of TAD
            if TAD_boundary_kawagutii(j)-dist_to_TAD > 0 & TAD_boundary_kawagutii(j)+dist_to_TAD < chromosome(end,1)

                left_avg_TPM = median(mov_avg_TPM(TAD_boundary_kawagutii(j)-dist_to_TAD:TAD_boundary_kawagutii(j))');
                right_avg_TPM = median(mov_avg_TPM(TAD_boundary_kawagutii(j):TAD_boundary_kawagutii(j)+dist_to_TAD)');

                if right_avg_TPM > left_avg_TPM
                    TAD_TPM_kawagutii = [TAD_TPM_kawagutii; -1.*(mov_avg_TPM(TAD_boundary_kawagutii(j)-dist_to_TAD:TAD_boundary_kawagutii(j)+dist_to_TAD)')];
                else
                    TAD_TPM_kawagutii = [TAD_TPM_kawagutii; mov_avg_TPM(TAD_boundary_kawagutii(j)-dist_to_TAD:TAD_boundary_kawagutii(j)+dist_to_TAD)'];
                end
                %elseif TAD_boundary_kawagutii(j)-dist_to_TAD < 0
                %TAD_TPM = [TAD_TPM; mov_avg_TPM(1:TAD_boundary_kawagutii(j)+dist_to_TAD)'];
                %elseif TAD_boundary_kawagutii(j)+dist_to_TAD > chromosome(end,1)
                %TAD_TPM = [TAD_TPM; mov_avg_TPM(TAD_boundary_kawagutii(j)-dist_to_TAD:chromosome(end,1))'];
            end
        end

        [coeff,score] = pca(chromosome(:,2:4));
        chromosome_PCA = score;

        P_axis_1_kawagutii=[P_axis_1_kawagutii; prctile(abs(score(:,1)),95)];
        P_axis_2_kawagutii=[P_axis_2_kawagutii; prctile(abs(score(:,2)),95)];
        P_axis_3_kawagutii=[P_axis_3_kawagutii; prctile(abs(score(:,3)),95)];

        for j = 1:1:size(TADs_kawagutii,1)
            TAD_PCA = chromosome_PCA(round(TADs_kawagutii(j,1)./HiC_resolution)+1:round(TADs_kawagutii(j,2)./HiC_resolution),:);
            %calculate moment of inertia tensor
            %Following: Chu, X., & Wang, J. (2023). Quantifying the large-scale chromosome structural dynamics during the mitosis-to-G1 phase transition of cell cycle. Open Biology, 13(11), 230175.
            Itensor_TAD = [[TAD_PCA(:,1)'*TAD_PCA(:,1) TAD_PCA(:,1)'*TAD_PCA(:,2) TAD_PCA(:,1)'*TAD_PCA(:,3)]; [TAD_PCA(:,2)'*TAD_PCA(:,1) TAD_PCA(:,2)'*TAD_PCA(:,2) TAD_PCA(:,2)'*TAD_PCA(:,3)]; [TAD_PCA(:,3)'*TAD_PCA(:,1) TAD_PCA(:,3)'*TAD_PCA(:,2) TAD_PCA(:,3)'*TAD_PCA(:,3)]];
            eigenvalues = eig(Itensor_TAD);
            asphericity_tensor_TADs_kawagutii = [asphericity_tensor_TADs_kawagutii; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];
        end
        Itensor_chromosome = [[chromosome_PCA(:,1)'*chromosome_PCA(:,1) chromosome_PCA(:,1)'*chromosome_PCA(:,2) chromosome_PCA(:,1)'*chromosome_PCA(:,3)]; [chromosome_PCA(:,2)'*chromosome_PCA(:,1) chromosome_PCA(:,2)'*chromosome_PCA(:,2) chromosome_PCA(:,2)'*chromosome_PCA(:,3)]; [chromosome_PCA(:,3)'*chromosome_PCA(:,1) chromosome_PCA(:,3)'*chromosome_PCA(:,2) chromosome_PCA(:,3)'*chromosome_PCA(:,3)]];
        eigenvalues = eig(Itensor_chromosome);
        asphericity_tensor_kawagutii = [asphericity_tensor_kawagutii; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];
        
        i
    else
        % File does not exist.
    end
end

%TPM versus distance from TAD boundary
figure
hold on
plot((-dist_to_TAD:1:dist_to_TAD)./1000,mean(TAD_TPM_kawagutii,1), 'LineWidth', 8, 'Color', [0.4660 0.6740 0.1880])
xline(0)
yline(0)
ylim([-8 6])
ax = gca;
ax.FontSize = 16;
xlabel('Distance from TAD Boundary [kbp]','FontSize', 24)
ylabel('TPM','FontSize', 24)

%TAD asphericity
figure
hold on
histogram(asphericity_tensor_TADs_kawagutii, 20, FaceColor = [0.4660 0.6740 0.1880]./3)  %a lighter colour
histogram(asphericity_tensor_kawagutii, FaceColor = [0.4660 0.6740 0.1880])
lgd=legend({'TADs Only','Entire Chromosome'});
legend boxoff
hold off
xlim([0 1])
ax = gca;
ax.FontSize = 16;
xlabel('Asphericity','FontSize', 24)
ylabel('Count','FontSize', 24)

% %correlate # of gene orientation switches with chromosome asphericity
% 
% %correlate # of TADs with chromosome asphericity
% figure
% hold on
% scatter(TAD_num_microadriaticum, asphericity_tensor_microadriaticum, 60, 'blue', 'filled')
% scatter(TAD_num_kawagutii, asphericity_tensor_kawagutii, 60, 'green', 'filled')
% lgd=legend({sprintf('Symbiodinium microadriaticum'),sprintf('Symbiodinium kawagutii')});
% legend boxoff
% hold off
% ax = gca;
% ax.FontSize = 16;
% %ylim([1 2])
% xlabel('Number of TADs','FontSize', 24)
% ylabel('Chromosome Asphericity','FontSize', 24)
% 
% %add 0 TAD chromosomes to scatter plots
% 
% figure
% hold on
% scatter(TAD_percent_coverage_microadriaticum, asphericity_tensor_microadriaticum, 25, 'blue', 'filled')
% scatter(TAD_percent_coverage_kawagutii, asphericity_tensor_kawagutii, 25, 'green', 'filled')
% lgd=legend({sprintf('Symbiodinium microadriaticum'),sprintf('Symbiodinium kawagutii')});
% legend boxoff
% hold off
% ax = gca;
% ax.FontSize = 16;
% %ylim([1 2])
% xlabel('TAD % Coverage','FontSize', 24)
% ylabel('Chromosome Asphericity','FontSize', 24)

