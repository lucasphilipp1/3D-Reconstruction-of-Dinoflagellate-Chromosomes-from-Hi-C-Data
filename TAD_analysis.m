%Fig 4 B in paper: histograms of TAD asphericity
%Fig S11 in paper: convergent transcription at TAD boundaries

%also writes chr #, TAD start, TAD end, TAD asphericity tensor, TAD asphericity PCA
%to a text file

clc
clear

num_chroms = 50;
dist_to_TAD = 45000;
bin_size = 1000;

asphericity_tensor_TADs_microadriaticum =[];
asphericity_tensor_chromosome_microadriaticum = [];

asphericity_tensor_TADs_kawagutii=[];
asphericity_tensor_chromosome_kawagutii = [];

asphericity_PCA_TADs_microadriaticum = [];
asphericity_PCA_chromosome_microadriaticum = [];

asphericity_PCA_TADs_kawagutii = [];
asphericity_PCA_chromosome_kawagutii = [];

asphericity_PCA_TADs_microadriaticum_write = {};
asphericity_tensor_TADs_microadriaticum_write = {};

asphericity_PCA_TADs_kawagutii_write = {};
asphericity_tensor_TADs_kawagutii_write = {};

%want to separately average convergent and divergent transcription at TAD boundaries
%don't want negative and positive TPM values to cancel out
TAD_TPM_left_positive_microadriaticum = [];
TAD_TPM_right_positive_microadriaticum = [];
TAD_TPM_left_positive_kawagutii = [];
TAD_TPM_right_positive_kawagutii = [];

P_axis_1_microadriaticum = [];
P_axis_2_microadriaticum = [];
P_axis_3_microadriaticum = [];

P_axis_1_kawagutii = [];
P_axis_2_kawagutii = [];
P_axis_3_kawagutii = [];

for i=1:num_chroms
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
            temp_TAD_boundary=[temp_TAD_boundary; str2double(extractBetween(TADfile{j},istab(1)+1,istab(2)-1)) str2double(extractBetween(TADfile{j},istab(2)+1,istab(3)-1))];
        end

        TADs_microadriaticum = [TADs_microadriaticum; temp_TAD_boundary];

        TPM_bp = zeros(chromosome(end,1),1);
        for j = 1:1:size(genes_start_end_TPM,1)
            TPM_bp(genes_start_end_TPM(j,1):genes_start_end_TPM(j,2)) = genes_start_end_TPM(j,3);
        end

        %sliding window
        mov_avg_TPM = movmedian(TPM_bp,bin_size);

        %left TAD boundary
        for j = 1:1:size(TADs_microadriaticum,1)
            %look dist_to_TAD bp upstream and downstream of TAD
            %ensure left TAD boundary is not too close to start or end of chromosome
            if TADs_microadriaticum(j,1)-dist_to_TAD > 0 & TADs_microadriaticum(j,1)+dist_to_TAD < chromosome(end,1)

                left_avg_TPM = median(mov_avg_TPM(TADs_microadriaticum(j,1)-dist_to_TAD:TADs_microadriaticum(j,1))');
                right_avg_TPM = median(mov_avg_TPM(TADs_microadriaticum(j,1):TADs_microadriaticum(j,1)+dist_to_TAD)');

                if right_avg_TPM > left_avg_TPM
                    TAD_TPM_right_positive_microadriaticum = [TAD_TPM_right_positive_microadriaticum; mov_avg_TPM(TADs_microadriaticum(j,1)-dist_to_TAD:TADs_microadriaticum(j,1)+dist_to_TAD)'];
                else
                    TAD_TPM_left_positive_microadriaticum = [TAD_TPM_left_positive_microadriaticum; mov_avg_TPM(TADs_microadriaticum(j,1)-dist_to_TAD:TADs_microadriaticum(j,1)+dist_to_TAD)'];
                end
            end
        end

        %right TAD boundary
        for j = 1:1:size(TADs_microadriaticum,1)
            %look dist_to_TAD bp upstream and downstream of TAD
            %ensure right TAD boundary is not too close to start or end of chromosome            
            if TADs_microadriaticum(j,2)-dist_to_TAD > 0 & TADs_microadriaticum(j,2)+dist_to_TAD < chromosome(end,1)

                left_avg_TPM = median(mov_avg_TPM(TADs_microadriaticum(j,2)-dist_to_TAD:TADs_microadriaticum(j,2))');
                right_avg_TPM = median(mov_avg_TPM(TADs_microadriaticum(j,2):TADs_microadriaticum(j,2)+dist_to_TAD)');

                if right_avg_TPM > left_avg_TPM
                    TAD_TPM_right_positive_microadriaticum = [TAD_TPM_right_positive_microadriaticum; mov_avg_TPM(TADs_microadriaticum(j,2)-dist_to_TAD:TADs_microadriaticum(j,2)+dist_to_TAD)'];
                else
                    TAD_TPM_left_positive_microadriaticum = [TAD_TPM_left_positive_microadriaticum; mov_avg_TPM(TADs_microadriaticum(j,2)-dist_to_TAD:TADs_microadriaticum(j,2)+dist_to_TAD)'];
                end
            end
        end

        [coeff,score] = pca(chromosome(:,2:4));
        chromosome_PCA = score;

        for j = 1:1:size(TADs_microadriaticum,1)
            TAD_PCA = chromosome_PCA(round(TADs_microadriaticum(j,1)./HiC_resolution)+1:round(TADs_microadriaticum(j,2)./HiC_resolution),:);
            %calculate moment of inertia tensor
            %Following: Huang, W., & Zaburdaev, V. (2019). The shape of pinned forced polymer loops. Soft Matter, 15(8), 1785-1792.

            TAD_PCA(:,1) = TAD_PCA(:,1) - mean(TAD_PCA(:,1));
            TAD_PCA(:,2) = TAD_PCA(:,2) - mean(TAD_PCA(:,2));
            TAD_PCA(:,3) = TAD_PCA(:,3) - mean(TAD_PCA(:,3));

            [coeff,score] = pca(TAD_PCA);
            TAD_PCA = score;

            xcm = sum(TAD_PCA(:,1))/(size(TAD_PCA,1));
            ycm = sum(TAD_PCA(:,2))/(size(TAD_PCA,1));
            zcm = sum(TAD_PCA(:,3))/(size(TAD_PCA,1));
            Gyration_tensor = (1/size(TAD_PCA,1)).*[[sum((TAD_PCA(:,1)-xcm).^2) sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,2)-ycm)) sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,3)-zcm))]; [sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,2)-ycm)) sum((TAD_PCA(:,2)-ycm).^2) sum((TAD_PCA(:,2)-ycm).*(TAD_PCA(:,3)-zcm))]; [sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,3)-zcm)) sum((TAD_PCA(:,2)-ycm).*(TAD_PCA(:,3)-zcm)) sum((TAD_PCA(:,3)-zcm).^2)]];
            eigenvalues = eig(Gyration_tensor);
            eigenvalues=sort(eigenvalues,'descend');
            asphericity_tensor_TADs_microadriaticum = [asphericity_tensor_TADs_microadriaticum; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];
            asphericity_tensor_TADs_microadriaticum_write{end+1,1} = sprintf('chr%i_pilon',i);
            asphericity_tensor_TADs_microadriaticum_write{end,2} = TADs_microadriaticum(j,1);
            asphericity_tensor_TADs_microadriaticum_write{end,3} = TADs_microadriaticum(j,2);
            asphericity_tensor_TADs_microadriaticum_write{end,4} = (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2);

            P_axis_1=prctile(TAD_PCA(:,1)-min(TAD_PCA(:,1)),95);
            P_axis_2=prctile(TAD_PCA(:,2)-min(TAD_PCA(:,2)),95);
            P_axis_3=prctile(TAD_PCA(:,3)-min(TAD_PCA(:,3)),95);

            asphericity_PCA_TADs_microadriaticum = [asphericity_PCA_TADs_microadriaticum; ((P_axis_1/P_axis_2)+(P_axis_1/P_axis_3)+(P_axis_2/P_axis_3))/3];
            asphericity_PCA_TADs_microadriaticum_write{end+1,1} = sprintf('chr%i_pilon',i);
            asphericity_PCA_TADs_microadriaticum_write{end,2} = TADs_microadriaticum(j,1);
            asphericity_PCA_TADs_microadriaticum_write{end,3} = TADs_microadriaticum(j,2);
            asphericity_PCA_TADs_microadriaticum_write{end,4} = ((P_axis_1/P_axis_2)+(P_axis_1/P_axis_3)+(P_axis_2/P_axis_3))/3;
        end
        xcm = sum(chromosome_PCA(:,1))/(size(chromosome_PCA,1));
        ycm = sum(chromosome_PCA(:,2))/(size(chromosome_PCA,1));
        zcm = sum(chromosome_PCA(:,3))/(size(chromosome_PCA,1));
        Gyration_tensor = (1/size(chromosome_PCA,1)).*[[sum((chromosome_PCA(:,1)-xcm).^2) sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,2)-ycm)) sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,3)-zcm))]; [sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,2)-ycm)) sum((chromosome_PCA(:,2)-ycm).^2) sum((chromosome_PCA(:,2)-ycm).*(chromosome_PCA(:,3)-zcm))]; [sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,3)-zcm)) sum((chromosome_PCA(:,2)-ycm).*(chromosome_PCA(:,3)-zcm)) sum((chromosome_PCA(:,3)-zcm).^2)]];
        eigenvalues = eig(Gyration_tensor);
        eigenvalues=sort(eigenvalues,'descend');
        asphericity_tensor_chromosome_microadriaticum = [asphericity_tensor_chromosome_microadriaticum; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];

        P_axis_1=prctile(chromosome_PCA(:,1)-min(chromosome_PCA(:,1)),95);
        P_axis_2=prctile(chromosome_PCA(:,2)-min(chromosome_PCA(:,2)),95);
        P_axis_3=prctile(chromosome_PCA(:,3)-min(chromosome_PCA(:,3)),95);

        asphericity_PCA_chromosome_microadriaticum = [asphericity_PCA_chromosome_microadriaticum; ((P_axis_1/P_axis_2)+(P_axis_1/P_axis_3)+(P_axis_2/P_axis_3))/3];

        i
    else
        % File does not exist.
    end
end

%TPM versus distance from left TAD boundary
figure
hold on
plot((-dist_to_TAD:1:dist_to_TAD)./1000,mean(TAD_TPM_right_positive_microadriaticum,1), 'LineWidth', 3, 'Color', [0 0.4470 0.7410])
plot((-dist_to_TAD:1:dist_to_TAD)./1000,mean(TAD_TPM_left_positive_microadriaticum,1), 'LineWidth', 8, 'Color', [0 0.4470 0.7410])
xline(0)
yline(0)
ylim([-10 10])
ax = gca;
ax.FontSize = 16;
xlabel('Distance from TAD Boundary [kbp]','FontSize', 24)
ylabel('TPM','FontSize', 24)

%TAD asphericity
figure
hold on
%histogram(all_asphericity_tensor(:,4,1:2,1:2), 20, FaceColor = [0.9290 0.6940 0.1250])
histogram(asphericity_tensor_TADs_microadriaticum, linspace(0,1,50), FaceColor = [0 0.4470 0.7410]./3) %a lighter colour
histogram(asphericity_tensor_chromosome_microadriaticum, linspace(0,1,50), FaceColor = [0 0.4470 0.7410])
lgd=legend({'Equilibrium Globule Regions','TADs Only','Entire Chromosome'});
legend boxoff
hold off
xlim([0 1])
ax = gca;
ax.FontSize = 16;
xlabel('Asphericity (tensor)','FontSize', 24)
ylabel('Count','FontSize', 24)

%TAD asphericity
figure
hold on
histogram(asphericity_PCA_TADs_microadriaticum, linspace(1,3,20), FaceColor = [0 0.4470 0.7410]./3) %a lighter colour
histogram(asphericity_PCA_chromosome_microadriaticum, FaceColor = [0 0.4470 0.7410])
lgd=legend({'TADs Only','Entire Chromosome'});
legend boxoff
hold off
ax = gca;
ax.FontSize = 16;
xlabel('Asphericity (PCA)','FontSize', 24)
ylabel('Count','FontSize', 24)
xlim([1 3])

for i=1:num_chroms
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
            temp_TAD_boundary=[temp_TAD_boundary; str2double(extractBetween(TADfile{j},istab(1)+1,istab(2)-1)) str2double(extractBetween(TADfile{j},istab(2)+1,istab(3)-1))];
        end

        TADs_kawagutii = [TADs_kawagutii; temp_TAD_boundary];

        TPM_bp = zeros(chromosome(end,1),1);
        for j = 1:1:size(genes_start_end_TPM,1)
            TPM_bp(genes_start_end_TPM(j,1):genes_start_end_TPM(j,2)) = genes_start_end_TPM(j,3);
        end

        %sliding window
        mov_avg_TPM = movmedian(TPM_bp,bin_size);

        %left TAD boundary
        for j = 1:1:size(TADs_kawagutii,1)
            %look dist_to_TAD bp upstream and downstream of TAD
            %ensure left TAD boundary is not too close to start or end of chromosome
            if TADs_kawagutii(j,1)-dist_to_TAD > 0 & TADs_kawagutii(j,1)+dist_to_TAD < chromosome(end,1)

                left_avg_TPM = median(mov_avg_TPM(TADs_kawagutii(j,1)-dist_to_TAD:TADs_kawagutii(j,1))');
                right_avg_TPM = median(mov_avg_TPM(TADs_kawagutii(j,1):TADs_kawagutii(j,1)+dist_to_TAD)');

                if right_avg_TPM > left_avg_TPM
                    TAD_TPM_right_positive_kawagutii = [TAD_TPM_right_positive_kawagutii; mov_avg_TPM(TADs_kawagutii(j,1)-dist_to_TAD:TADs_kawagutii(j,1)+dist_to_TAD)'];
                else
                    TAD_TPM_left_positive_kawagutii = [TAD_TPM_left_positive_kawagutii; mov_avg_TPM(TADs_kawagutii(j,1)-dist_to_TAD:TADs_kawagutii(j,1)+dist_to_TAD)'];
                end
            end
        end

        %right TAD boundary
        for j = 1:1:size(TADs_kawagutii,1)
            %look dist_to_TAD bp upstream and downstream of TAD
            %ensure left TAD boundary is not too close to start or end of chromosome
            if TADs_kawagutii(j,2)-dist_to_TAD > 0 & TADs_kawagutii(j,2)+dist_to_TAD < chromosome(end,1)

                left_avg_TPM = median(mov_avg_TPM(TADs_kawagutii(j,2)-dist_to_TAD:TADs_kawagutii(j,2))');
                right_avg_TPM = median(mov_avg_TPM(TADs_kawagutii(j,2):TADs_kawagutii(j,2)+dist_to_TAD)');

                if right_avg_TPM > left_avg_TPM
                    TAD_TPM_right_positive_kawagutii = [TAD_TPM_right_positive_kawagutii; mov_avg_TPM(TADs_kawagutii(j,2)-dist_to_TAD:TADs_kawagutii(j,2)+dist_to_TAD)'];
                else
                    TAD_TPM_left_positive_kawagutii = [TAD_TPM_left_positive_kawagutii; mov_avg_TPM(TADs_kawagutii(j,2)-dist_to_TAD:TADs_kawagutii(j,2)+dist_to_TAD)'];
                end
            end
        end

        [coeff,score] = pca(chromosome(:,2:4));
        chromosome_PCA = score;

        for j = 1:1:size(TADs_kawagutii,1)
            TAD_PCA = chromosome_PCA(round(TADs_kawagutii(j,1)./HiC_resolution)+1:round(TADs_kawagutii(j,2)./HiC_resolution),:);
            %calculate moment of inertia tensor
            %Following: Huang, W., & Zaburdaev, V. (2019). The shape of pinned forced polymer loops. Soft Matter, 15(8), 1785-1792.
            
            TAD_PCA(:,1) = TAD_PCA(:,1) - mean(TAD_PCA(:,1));
            TAD_PCA(:,2) = TAD_PCA(:,2) - mean(TAD_PCA(:,2));
            TAD_PCA(:,3) = TAD_PCA(:,3) - mean(TAD_PCA(:,3));

            [coeff,score] = pca(TAD_PCA);
            TAD_PCA = score;

            xcm = sum(TAD_PCA(:,1))/(size(TAD_PCA,1));
            ycm = sum(TAD_PCA(:,2))/(size(TAD_PCA,1));
            zcm = sum(TAD_PCA(:,3))/(size(TAD_PCA,1));
            Gyration_tensor = (1/size(TAD_PCA,1)).*[[sum((TAD_PCA(:,1)-xcm).^2) sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,2)-ycm)) sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,3)-zcm))]; [sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,2)-ycm)) sum((TAD_PCA(:,2)-ycm).^2) sum((TAD_PCA(:,2)-ycm).*(TAD_PCA(:,3)-zcm))]; [sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,3)-zcm)) sum((TAD_PCA(:,2)-ycm).*(TAD_PCA(:,3)-zcm)) sum((TAD_PCA(:,3)-zcm).^2)]];
            eigenvalues = eig(Gyration_tensor);
            eigenvalues=sort(eigenvalues,'descend');
            asphericity_tensor_TADs_kawagutii = [asphericity_tensor_TADs_kawagutii; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];
            asphericity_tensor_TADs_kawagutii_write{end+1,1} = sprintf('HiC_scaffold_%i',i);
            asphericity_tensor_TADs_kawagutii_write{end,2} = TADs_kawagutii(j,1);
            asphericity_tensor_TADs_kawagutii_write{end,3} = TADs_kawagutii(j,2);
            asphericity_tensor_TADs_kawagutii_write{end,4} = (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2);

            P_axis_1=prctile(TAD_PCA(:,1)-min(TAD_PCA(:,1)),95);
            P_axis_2=prctile(TAD_PCA(:,2)-min(TAD_PCA(:,2)),95);
            P_axis_3=prctile(TAD_PCA(:,3)-min(TAD_PCA(:,3)),95);

            asphericity_PCA_TADs_kawagutii = [asphericity_PCA_TADs_kawagutii; ((P_axis_1/P_axis_2)+(P_axis_1/P_axis_3)+(P_axis_2/P_axis_3))/3];
            asphericity_PCA_TADs_kawagutii_write{end+1,1} = sprintf('HiC_scaffold_%i',i);
            asphericity_PCA_TADs_kawagutii_write{end,2} = TADs_kawagutii(j,1);
            asphericity_PCA_TADs_kawagutii_write{end,3} = TADs_kawagutii(j,2);
            asphericity_PCA_TADs_kawagutii_write{end,4} = ((P_axis_1/P_axis_2)+(P_axis_1/P_axis_3)+(P_axis_2/P_axis_3))/3;
        end
        xcm = sum(chromosome_PCA(:,1))/(size(chromosome_PCA,1));
        ycm = sum(chromosome_PCA(:,2))/(size(chromosome_PCA,1));
        zcm = sum(chromosome_PCA(:,3))/(size(chromosome_PCA,1));
        Gyration_tensor = (1/size(chromosome_PCA,1)).*[[sum((chromosome_PCA(:,1)-xcm).^2) sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,2)-ycm)) sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,3)-zcm))]; [sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,2)-ycm)) sum((chromosome_PCA(:,2)-ycm).^2) sum((chromosome_PCA(:,2)-ycm).*(chromosome_PCA(:,3)-zcm))]; [sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,3)-zcm)) sum((chromosome_PCA(:,2)-ycm).*(chromosome_PCA(:,3)-zcm)) sum((chromosome_PCA(:,3)-zcm).^2)]];
        eigenvalues = eig(Gyration_tensor);
        eigenvalues=sort(eigenvalues,'descend');
        asphericity_tensor_chromosome_kawagutii = [asphericity_tensor_chromosome_kawagutii; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];

        P_axis_1=prctile(chromosome_PCA(:,1)-min(chromosome_PCA(:,1)),95);
        P_axis_2=prctile(chromosome_PCA(:,2)-min(chromosome_PCA(:,2)),95);
        P_axis_3=prctile(chromosome_PCA(:,3)-min(chromosome_PCA(:,3)),95);

        asphericity_PCA_chromosome_kawagutii = [asphericity_PCA_chromosome_kawagutii; ((P_axis_1/P_axis_2)+(P_axis_1/P_axis_3)+(P_axis_2/P_axis_3))/3];

        i
    else
        % File does not exist.
    end
end

%TPM versus distance from TAD boundary
figure
hold on
plot((-dist_to_TAD:1:dist_to_TAD)./1000,mean(TAD_TPM_right_positive_kawagutii,1), 'LineWidth', 3, 'Color', [0.4660 0.6740 0.1880])
plot((-dist_to_TAD:1:dist_to_TAD)./1000,mean(TAD_TPM_left_positive_kawagutii,1), 'LineWidth', 8, 'Color', [0.4660 0.6740 0.1880])
xline(0)
yline(0)
ylim([-10 10])
ax = gca;
ax.FontSize = 16;
xlabel('Distance from TAD Boundary [kbp]','FontSize', 24)
ylabel('TPM','FontSize', 24)

%TAD asphericity
figure
hold on 
%histogram(all_asphericity_tensor(:,4,1:2,1:2), 20, FaceColor = [0.9290 0.6940 0.1250])
histogram(asphericity_tensor_TADs_kawagutii, linspace(0,1,50), FaceColor = [0.4660 0.6740 0.1880]./3)  %a lighter colour
histogram(asphericity_tensor_chromosome_kawagutii, linspace(0,1,50), FaceColor = [0.4660 0.6740 0.1880])
lgd=legend({'Equilibrium Globule Regions','TADs Only','Entire Chromosome'});
legend boxoff
hold off
xlim([0 1])
ax = gca;
ax.FontSize = 16;
xlabel('Asphericity (tensor)','FontSize', 24)
ylabel('Count','FontSize', 24)

%TAD asphericity
figure
hold on
histogram(asphericity_PCA_TADs_kawagutii, linspace(1,3,20), FaceColor = [0.4660 0.6740 0.1880]./3)  %a lighter colour
histogram(asphericity_PCA_chromosome_kawagutii, FaceColor = [0.4660 0.6740 0.1880])
lgd=legend({'TADs Only','Entire Chromosome'});
legend boxoff
hold off
ax = gca;
ax.FontSize = 16;
xlabel('Asphericity (PCA)','FontSize', 24)
ylabel('Count','FontSize', 24)
xlim([1 3])

writecell(asphericity_PCA_TADs_microadriaticum_write,'symbiodinium_microadriaticum_allchroms_TADs_asphericity_PCA.txt','Delimiter','\t')
writecell(asphericity_PCA_TADs_kawagutii_write,'symbiodinium_kawagutii_allchroms_TADs_asphericity_PCA.txt','Delimiter','\t')

writecell(asphericity_tensor_TADs_microadriaticum_write,'symbiodinium_microadriaticum_allchroms_TADs_asphericity_tensor.txt','Delimiter','\t')
writecell(asphericity_tensor_TADs_kawagutii_write,'symbiodinium_kawagutii_allchroms_TADs_asphericity_tensor.txt','Delimiter','\t')