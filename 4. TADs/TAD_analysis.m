%Fig 3 D in paper: histograms of TAD asphericity
%Fig S10 in paper: convergent transcription at TAD boundaries

num_chroms = 50;
dist_to_TAD = 45000;
bin_size = 1000;

[TAD_TPM_right_positive_microadriaticum, TAD_TPM_left_positive_microadriaticum asphericity_tensor_TADs_microadriaticum asphericity_tensor_chromosome_microadriaticum asphericity_PCA_TADs_microadriaticum asphericity_PCA_chromosome_microadriaticum] = gene_expression_vs_dist_to_TAD_boundary('symbiodinium_microadriaticum_chr%i_3D.xyz', 's_microadriaticum_chr%i_pilon_strandsum_TPM.bed', 'symbiodinium_microadriaticum_chr%i.txt.bed', num_chroms, dist_to_TAD, bin_size);

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

[TAD_TPM_right_positive_kawagutii TAD_TPM_left_positive_kawagutii asphericity_tensor_TADs_kawagutii asphericity_tensor_chromosome_kawagutii asphericity_PCA_TADs_kawagutii asphericity_PCA_chromosome_kawagutii] = gene_expression_vs_dist_to_TAD_boundary('s_kawagutii_V3_HiC_scaffold_%i.xyz', 's_kawagutii_chr%i_strandsum_TPM.bed', 's_kawagutii_V3_HiC_scaffold_%i_for_CSynth.txt.bed', num_chroms, dist_to_TAD, bin_size);

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

%250mon/disc*5000bp/mon = 16 discs
%150mon/disc*5000bp/mon = 27 discs
%75mon/disc*5000bp/mon = 54 discs

asphericity_tensor_disc = [];
asphericity_PCA_disc = [];

asphericity_tensor_CLC_16_discs = [];
asphericity_tensor_CLC_27_discs = [];
asphericity_tensor_CLC_54_discs = [];

for k = 1:1:3
    for i = 1:num_chroms
        if k == 1
            chromosome = importdata(sprintf('cholesteric_monomer_locations_16_discs_%i.txt',i));
            chromosome = [(1:1:size(chromosome,1))' chromosome];
            ndiscs = 16;
        elseif k == 2
            chromosome = importdata(sprintf('cholesteric_monomer_locations_27_discs_%i.txt',i));
            chromosome = [(1:1:size(chromosome,1))' chromosome];
            ndiscs = 27;
        elseif k == 3
            chromosome = importdata(sprintf('cholesteric_monomer_locations_54_discs_%i.txt',i));
            chromosome = [(1:1:size(chromosome,1))' chromosome];
            ndiscs = 54;
        end

        chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
        chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));
        chromosome(:,4) = chromosome(:,4) - mean(chromosome(:,4));

        %break up into individual disks
        for i = 1:1:ndiscs
            [coeff,score] = pca(chromosome(1+(i-1)*ndiscs:i*ndiscs,2:4));
            disc_PCA = score;

            disc_PCA(:,1) = disc_PCA(:,1) - mean(disc_PCA(:,1));
            disc_PCA(:,2) = disc_PCA(:,2) - mean(disc_PCA(:,2));
            disc_PCA(:,3) = disc_PCA(:,3) - mean(disc_PCA(:,3));

            xcm = sum(disc_PCA(:,1))/(size(disc_PCA,1));
            ycm = sum(disc_PCA(:,2))/(size(disc_PCA,1));
            zcm = sum(disc_PCA(:,3))/(size(disc_PCA,1));
            Gyration_tensor = (1/size(disc_PCA,1)).*[[sum((disc_PCA(:,1)-xcm).^2) sum((disc_PCA(:,1)-xcm).*(disc_PCA(:,2)-ycm)) sum((disc_PCA(:,1)-xcm).*(disc_PCA(:,3)-zcm))]; [sum((disc_PCA(:,1)-xcm).*(disc_PCA(:,2)-ycm)) sum((disc_PCA(:,2)-ycm).^2) sum((disc_PCA(:,2)-ycm).*(disc_PCA(:,3)-zcm))]; [sum((disc_PCA(:,1)-xcm).*(disc_PCA(:,3)-zcm)) sum((disc_PCA(:,2)-ycm).*(disc_PCA(:,3)-zcm)) sum((disc_PCA(:,3)-zcm).^2)]];
            eigenvalues = eig(Gyration_tensor);
            eigenvalues=sort(eigenvalues,'descend');
            asphericity_tensor_disc = [asphericity_tensor_disc; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];

            %average ratio of PCs
            P_axis_1=prctile(disc_PCA(:,1)-min(disc_PCA(:,1)),95);
            P_axis_2=prctile(disc_PCA(:,2)-min(disc_PCA(:,2)),95);
            P_axis_3=prctile(disc_PCA(:,3)-min(disc_PCA(:,3)),95);

            asphericity_PCA_disc = [asphericity_PCA_disc; ((P_axis_1/P_axis_2)+(P_axis_1/P_axis_3)+(P_axis_2/P_axis_3))/3];

        end
        if k == 1
            asphericity_tensor_CLC_16_discs = [asphericity_tensor_CLC_16_discs; asphericity_tensor_disc];
        elseif k == 2
            asphericity_tensor_CLC_27_discs = [asphericity_tensor_CLC_27_discs; asphericity_tensor_disc];
        elseif k == 3
            asphericity_tensor_CLC_54_discs = [asphericity_tensor_CLC_54_discs; asphericity_tensor_disc];
        end
    end
end

%TAD asphericity
edges = linspace(0, 1, 21);

[N1,e1]=histcounts(asphericity_tensor_TADs_kawagutii, edges);
[N2,e2]=histcounts(asphericity_tensor_TADs_microadriaticum, edges);

[N3,e3]=histcounts(asphericity_tensor_CLC_16_discs, edges);
[N4,e4]=histcounts(asphericity_tensor_CLC_27_discs, edges);
[N5,e5]=histcounts(asphericity_tensor_CLC_54_discs, edges);

e1 = e1(2:end) - (e1(2)-e1(1))/2;
e2 = e2(2:end) - (e2(2)-e2(1))/2;
e3 = e3(2:end) - (e3(2)-e3(1))/2;
e4 = e4(2:end) - (e4(2)-e4(1))/2;
e5 = e5(2:end) - (e5(2)-e5(1))/2;

figure
hold on
plot(e1,N1, '--', Color = [0.4660 0.6740 0.1880], LineWidth=3)
plot(e2,N2, '--', Color = [0 0.4470 0.7410], LineWidth=3)
plot(e3,N3./400, Color = [0.2 0.2 0.2], LineWidth=3) %reduce sample size for visualization purposes
plot(e4,N4./400, Color = [0.4 0.4 0.4], LineWidth=3)
plot(e5,N5./400, Color = [0.6 0.6 0.6], LineWidth=3)

lgd=legend({'TADs {\it S. kawagutii}','TADs {\it S. microadriaticum}','CLC 16 Discs', 'CLC 27 Discs', 'CLC 54 Discs'});
lgd.FontSize = 20;
legend boxoff
hold off
xlim([0 1])
ylim([0 100])
ax = gca;
ax.FontSize = 16;
xlabel('Asphericity','FontSize', 24)
ylabel('Count','FontSize', 24)

function [TAD_TPM_right_positive, TAD_TPM_left_positive asphericity_tensor_TADs asphericity_tensor_chromosome asphericity_PCA_TADs asphericity_PCA_chromosome] = gene_expression_vs_dist_to_TAD_boundary(chromosome_file, RNAseq_file, TAD_file, num_chroms, dist_to_TAD, bin_size)

asphericity_tensor_TADs=[];
asphericity_tensor_chromosome = [];

asphericity_PCA_TADs= [];
asphericity_PCA_chromosome = [];

%want to separately average convergent and divergent transcription at TAD boundaries
%don't want negative and positive TPM values to cancel out
TAD_TPM_left_positive = [];
TAD_TPM_right_positive = [];

for i=1:num_chroms
    TADs_collected = [];
    temp_TAD_boundary = [];

    if isfile(sprintf(TAD_file,i))
        % File exists.
        chromosome = importdata(sprintf(chromosome_file,i));
        RNAseq = importdata(sprintf(RNAseq_file,i));
        TAD = importdata(sprintf(TAD_file,i));

        chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
        chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));
        chromosome(:,4) = chromosome(:,4) - mean(chromosome(:,4));

        HiC_resolution = chromosome(2,1)-chromosome(1,1);

        genes_start_end_TPM = [RNAseq.data(:,1) RNAseq.data(:,2) RNAseq.data(:,8)];

        for j = 1:1:size(TAD,1)
            istab=strfind(TAD{j},char(9));
            temp_TAD_boundary=[temp_TAD_boundary; str2double(extractBetween(TAD{j},istab(1)+1,istab(2)-1)) str2double(extractBetween(TAD{j},istab(2)+1,istab(3)-1))];
        end

        TADs_collected = [TADs_collected; temp_TAD_boundary];

        TPM_bp = zeros(chromosome(end,1),1);
        for j = 1:1:size(genes_start_end_TPM,1)
            TPM_bp(genes_start_end_TPM(j,1):genes_start_end_TPM(j,2)) = genes_start_end_TPM(j,3);
        end

        %sliding window
        mov_avg_TPM = movmedian(TPM_bp,bin_size);


        %left TAD boundary
        for j = 1:1:size(TADs_collected,1)
            %look dist_to_TAD bp upstream and downstream of TAD
            %ensure left TAD boundary is not too close to start or end of chromosome
            if TADs_collected(j,1)-dist_to_TAD > 0 & TADs_collected(j,1)+dist_to_TAD < chromosome(end,1)

                left_avg_TPM = median(mov_avg_TPM(TADs_collected(j,1)-dist_to_TAD:TADs_collected(j,1))');
                right_avg_TPM = median(mov_avg_TPM(TADs_collected(j,1):TADs_collected(j,1)+dist_to_TAD)');

                %want to separately average convergent and divergent transcription at TAD boundaries
                %don't want negative and positive TPM values to cancel out

                if right_avg_TPM > left_avg_TPM
                    TAD_TPM_right_positive = [TAD_TPM_right_positive; mov_avg_TPM(TADs_collected(j,1)-dist_to_TAD:TADs_collected(j,1)+dist_to_TAD)'];
                else
                    TAD_TPM_left_positive = [TAD_TPM_left_positive; mov_avg_TPM(TADs_collected(j,1)-dist_to_TAD:TADs_collected(j,1)+dist_to_TAD)'];
                end
            end
        end

        %right TAD boundary
        for j = 1:1:size(TADs_collected,1)
            %look dist_to_TAD bp upstream and downstream of TAD
            %ensure left TAD boundary is not too close to start or end of chromosome
            if TADs_collected(j,2)-dist_to_TAD > 0 & TADs_collected(j,2)+dist_to_TAD < chromosome(end,1)

                left_avg_TPM = median(mov_avg_TPM(TADs_collected(j,2)-dist_to_TAD:TADs_collected(j,2))');
                right_avg_TPM = median(mov_avg_TPM(TADs_collected(j,2):TADs_collected(j,2)+dist_to_TAD)');

                %want to separately average convergent and divergent transcription at TAD boundaries
                %don't want negative and positive TPM values to cancel out

                if right_avg_TPM > left_avg_TPM
                    TAD_TPM_right_positive = [TAD_TPM_right_positive; mov_avg_TPM(TADs_collected(j,2)-dist_to_TAD:TADs_collected(j,2)+dist_to_TAD)'];
                else
                    TAD_TPM_left_positive = [TAD_TPM_left_positive; mov_avg_TPM(TADs_collected(j,2)-dist_to_TAD:TADs_collected(j,2)+dist_to_TAD)'];
                end
            end
        end

        [coeff,score] = pca(chromosome(:,2:4));
        chromosome_PCA = score;

        for j = 1:1:size(TADs_collected,1)
            TAD_PCA = chromosome_PCA(round(TADs_collected(j,1)./HiC_resolution)+1:round(TADs_collected(j,2)./HiC_resolution),:);
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

            R_g = (sum(TAD_PCA(:,1).^2+TAD_PCA(:,2).^2+TAD_PCA(:,3).^2)/size(TAD_PCA,1))^0.5; %radius of gyration
            [k v] = boundary(TAD_PCA(:,1),TAD_PCA(:,2),TAD_PCA(:,3));
            density = size(TAD_PCA,1)/v; %number of monomers/volume = DNA density

            Gyration_tensor = (1/size(TAD_PCA,1)).*[[sum((TAD_PCA(:,1)-xcm).^2) sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,2)-ycm)) sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,3)-zcm))]; [sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,2)-ycm)) sum((TAD_PCA(:,2)-ycm).^2) sum((TAD_PCA(:,2)-ycm).*(TAD_PCA(:,3)-zcm))]; [sum((TAD_PCA(:,1)-xcm).*(TAD_PCA(:,3)-zcm)) sum((TAD_PCA(:,2)-ycm).*(TAD_PCA(:,3)-zcm)) sum((TAD_PCA(:,3)-zcm).^2)]];
            eigenvalues = eig(Gyration_tensor);
            eigenvalues=sort(eigenvalues,'descend');
            asphericity_tensor_TADs = [asphericity_tensor_TADs; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];
            
            P_axis_1=prctile(TAD_PCA(:,1)-min(TAD_PCA(:,1)),95);
            P_axis_2=prctile(TAD_PCA(:,2)-min(TAD_PCA(:,2)),95);
            P_axis_3=prctile(TAD_PCA(:,3)-min(TAD_PCA(:,3)),95);

            asphericity_PCA_TADs = [asphericity_PCA_TADs; ((P_axis_1/P_axis_2)+(P_axis_1/P_axis_3)+(P_axis_2/P_axis_3))/3];
        end
        xcm = sum(chromosome_PCA(:,1))/(size(chromosome_PCA,1));
        ycm = sum(chromosome_PCA(:,2))/(size(chromosome_PCA,1));
        zcm = sum(chromosome_PCA(:,3))/(size(chromosome_PCA,1));
        Gyration_tensor = (1/size(chromosome_PCA,1)).*[[sum((chromosome_PCA(:,1)-xcm).^2) sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,2)-ycm)) sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,3)-zcm))]; [sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,2)-ycm)) sum((chromosome_PCA(:,2)-ycm).^2) sum((chromosome_PCA(:,2)-ycm).*(chromosome_PCA(:,3)-zcm))]; [sum((chromosome_PCA(:,1)-xcm).*(chromosome_PCA(:,3)-zcm)) sum((chromosome_PCA(:,2)-ycm).*(chromosome_PCA(:,3)-zcm)) sum((chromosome_PCA(:,3)-zcm).^2)]];
        eigenvalues = eig(Gyration_tensor);
        eigenvalues=sort(eigenvalues,'descend');
        asphericity_tensor_chromosome = [asphericity_tensor_chromosome; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];

        P_axis_1=prctile(chromosome_PCA(:,1)-min(chromosome_PCA(:,1)),95);
        P_axis_2=prctile(chromosome_PCA(:,2)-min(chromosome_PCA(:,2)),95);
        P_axis_3=prctile(chromosome_PCA(:,3)-min(chromosome_PCA(:,3)),95);

        asphericity_PCA_chromosome = [asphericity_PCA_chromosome; ((P_axis_1/P_axis_2)+(P_axis_1/P_axis_3)+(P_axis_2/P_axis_3))/3];

        i
    else
        % File does not exist.
    end
end
end
