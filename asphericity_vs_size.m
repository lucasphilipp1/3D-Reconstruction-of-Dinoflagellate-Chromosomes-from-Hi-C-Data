%Fig S13 in the paper
%quality control of asphericity metrics
%calculate asphericity for randomly selected regions of various sizes

clc
clear

num_chroms = 50; %don't increase past 50, smaller chromosomes won't be long enough for large region sizes

%define region sizes
region_sizes = linspace(10^4.5,0.75*10^7,250); %in units bp

N_sample_per_region_size = 10;
HiC_resolution = 5000; %bp per monomer

all_asphericity_tensor = zeros(size(region_sizes,2), 6, num_chroms, N_sample_per_region_size);
all_asphericity_PCA = zeros(size(region_sizes,2), 6, num_chroms, N_sample_per_region_size);

for l = 1:1:size(region_sizes,2)
    for k = 1:1:6
        for i = 1:num_chroms
            TADfile = [];
            TADs = [];
            temp_TAD_boundary = [];
            asphericity_tensor_rand_region = [];
            asphericity_PCA_rand_region = [];
            if k == 1
                chromosome = importdata(sprintf('cholesteric_monomer_locations_16_discs_%i.txt',i));
                chromosome = [(1:1:size(chromosome,1))'*HiC_resolution chromosome];
            elseif k == 2
                chromosome = importdata(sprintf('cholesteric_monomer_locations_27_discs_%i.txt',i));
                chromosome = [(1:1:size(chromosome,1))'*HiC_resolution chromosome];
            elseif k == 3
                chromosome = importdata(sprintf('cholesteric_monomer_locations_54_discs_%i.txt',i));
                chromosome = [(1:1:size(chromosome,1))'*HiC_resolution chromosome];
            elseif k == 4
                equilibrium_struct=tdfread(sprintf('equilibrium%i.dat',i));
                equilibrium_names = fieldnames(equilibrium_struct);
                steps_equil = size(getfield(equilibrium_struct,equilibrium_names{1}),1);
                equilibrium = zeros(steps_equil,3);
                chromosome = [(1:1:steps_equil)'*HiC_resolution getfield(equilibrium_struct,equilibrium_names{1})];
            elseif k == 5
                chromosome = importdata(sprintf('symbiodinium_microadriaticum_chr%i_3D.xyz',i));
                if isfile(sprintf('symbiodinium_microadriaticum_chr%i.txt.bed',i))
                    % File exists.
                    TADfile = importdata(sprintf('symbiodinium_microadriaticum_chr%i.txt.bed',i));
                else
                    % File does not exist.
                end
            elseif k == 6
                chromosome = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i.xyz',i));
                if isfile(sprintf('s_kawagutii_V3_HiC_scaffold_%i_for_CSynth.txt.bed',i))
                    % File exists.
                    TADfile = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i_for_CSynth.txt.bed',i));
                else
                    % File does not exist.
                end
            end

            chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
            chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));
            chromosome(:,4) = chromosome(:,4) - mean(chromosome(:,4));

            [coeff,score] = pca(chromosome(:,2:4));
            chromosome_PCA = score;

            starts = [];
            ends = [];

            starts = round(unifrnd(1,size(chromosome,1)-(round(region_sizes(l)/HiC_resolution)), [N_sample_per_region_size 1])); %in units monommer index
            ends = starts + round(region_sizes(l)/HiC_resolution); %in units monommer index

            %calculate asphericity of randomly selected regions
            for j = 1:1:size(starts,1)
                rand_region_PCA = chromosome_PCA(starts(j):ends(j),:);

                rand_region_PCA(:,1) = rand_region_PCA(:,1) - mean(rand_region_PCA(:,1));
                rand_region_PCA(:,2) = rand_region_PCA(:,2) - mean(rand_region_PCA(:,2));
                rand_region_PCA(:,3) = rand_region_PCA(:,3) - mean(rand_region_PCA(:,3));

                [coeff,score] = pca(rand_region_PCA);
                rand_region_PCA = score;

                xcm = sum(rand_region_PCA(:,1))/(size(rand_region_PCA,1));
                ycm = sum(rand_region_PCA(:,2))/(size(rand_region_PCA,1));
                zcm = sum(rand_region_PCA(:,3))/(size(rand_region_PCA,1));
                Gyration_tensor = (1/size(rand_region_PCA,1)).*[[sum((rand_region_PCA(:,1)-xcm).^2) sum((rand_region_PCA(:,1)-xcm).*(rand_region_PCA(:,2)-ycm)) sum((rand_region_PCA(:,1)-xcm).*(rand_region_PCA(:,3)-zcm))]; [sum((rand_region_PCA(:,1)-xcm).*(rand_region_PCA(:,2)-ycm)) sum((rand_region_PCA(:,2)-ycm).^2) sum((rand_region_PCA(:,2)-ycm).*(rand_region_PCA(:,3)-zcm))]; [sum((rand_region_PCA(:,1)-xcm).*(rand_region_PCA(:,3)-zcm)) sum((rand_region_PCA(:,2)-ycm).*(rand_region_PCA(:,3)-zcm)) sum((rand_region_PCA(:,3)-zcm).^2)]];
                eigenvalues = eig(Gyration_tensor);
                eigenvalues=sort(eigenvalues,'descend');
                asphericity_tensor_rand_region = [asphericity_tensor_rand_region; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];

                %average ratio of PCs
                P_axis_1=prctile(rand_region_PCA(:,1)-min(rand_region_PCA(:,1)),95);
                P_axis_2=prctile(rand_region_PCA(:,2)-min(rand_region_PCA(:,2)),95);
                P_axis_3=prctile(rand_region_PCA(:,3)-min(rand_region_PCA(:,3)),95);

                asphericity_PCA_rand_region = [asphericity_PCA_rand_region; ((P_axis_1/P_axis_2)+(P_axis_1/P_axis_3)+(P_axis_2/P_axis_3))/3];

            end

            all_asphericity_tensor(l,k,i,:) =  asphericity_tensor_rand_region;
            all_asphericity_PCA(l,k,i,:) =  asphericity_PCA_rand_region;

        end
        k
    end
    l
end

%collect all TAD sizes 
TAD_size_microadriaticum = [];
TAD_size_kawagutii = [];

for i=1:num_chroms
    TADs_microadriaticum = [];
    temp_TAD_boundary = [];

    if isfile(sprintf('symbiodinium_microadriaticum_chr%i.txt.bed',i))
        % File exists.
        TADfile = importdata(sprintf('symbiodinium_microadriaticum_chr%i.txt.bed',i));

        for j = 1:1:size(TADfile,1)
            istab=strfind(TADfile{j},char(9));
            temp_TAD_boundary=[temp_TAD_boundary; str2double(extractBetween(TADfile{j},istab(1)+1,istab(2)-1)) str2double(extractBetween(TADfile{j},istab(2)+1,istab(3)-1))];
        end

        TADs_microadriaticum = [TADs_microadriaticum; temp_TAD_boundary];

        for j = 1:1:size(TADs_microadriaticum,1)
            TAD_size_microadriaticum = [TAD_size_microadriaticum; TADs_microadriaticum(j,2)-TADs_microadriaticum(j,1)];
        end
        i
    else
        % File does not exist.
    end
end

for i=1:num_chroms
    TADs_kawagutii = [];
    temp_TAD_boundary = [];

    if isfile(sprintf('s_kawagutii_V3_HiC_scaffold_%i_for_CSynth.txt.bed',i))
        % File exists.
        TADfile = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i_for_CSynth.txt.bed',i));

        for j = 1:1:size(TADfile,1)
            istab=strfind(TADfile{j},char(9));
            temp_TAD_boundary=[temp_TAD_boundary; str2double(extractBetween(TADfile{j},istab(1)+1,istab(2)-1)) str2double(extractBetween(TADfile{j},istab(2)+1,istab(3)-1))];
        end

        TADs_kawagutii = [TADs_kawagutii; temp_TAD_boundary];

        for j = 1:1:size(TADs_kawagutii,1)
            TAD_size_kawagutii = [TAD_size_kawagutii; TADs_kawagutii(j,2)-TADs_kawagutii(j,1)];
        end
        i
    else
        % File does not exist.
    end
end

all_asphericity_tensor_mean = mean(all_asphericity_tensor, [3 4], 'omitnan');
all_asphericity_PCA_mean = mean(all_asphericity_PCA, [3 4], 'omitnan');

figure
hold on
plot(region_sizes, all_asphericity_tensor_mean(:,1), 'LineWidth',2, 'Color',[0 0 0])
plot(region_sizes, all_asphericity_tensor_mean(:,2), 'LineWidth',2, 'Color',[0.4 0.4 0.4])
plot(region_sizes, all_asphericity_tensor_mean(:,3), 'LineWidth',2, 'Color',[0.8 0.8 0.8])
plot(region_sizes, all_asphericity_tensor_mean(:,4), 'LineWidth',2, 'Color',[0.9290 0.6940 0.1250])
plot(region_sizes, all_asphericity_tensor_mean(:,5), 'LineWidth',2, 'Color', [0 0.4470 0.7410])
plot(region_sizes, all_asphericity_tensor_mean(:,6), 'LineWidth',2, 'Color', [0.4660 0.6740 0.1880])
xline((mean(asphericity_PCA_TADs_microadriaticum(:,1))+mean(asphericity_PCA_TADs_kawagutii(:,1)))/2,'--',{'Mean TAD Size'});
xlabel('Region Size [bp]', 'fontsize', 24)
ylabel('Asphericity (tensor)', 'fontsize', 24)
ax = gca;
ax.FontSize = 14;
legend('random regions CLC (16 discs)','random regions CLC (27 discs)','random regions CLC (54 discs)','random regions equilibrium globule', 'random regions S. microadriaticum', 'random regions S. kawagutii')
%legend('random regions CLC (16 discs)','random regions CLC (27 discs)','random regions CLC (54 discs)','random regions equilibrium globule', 'random regions S. microadriaticum', 'random regions S. kawagutii', 'whole chromosome S. microadriaticum', 'whole chromosome S. kawagutii', 'TADs S. microadriaticum', 'TADs S. kawagutii')
ylim([0 1])
legend boxoff
hold off

figure
hold on
plot(region_sizes, all_asphericity_PCA_mean(:,1), 'LineWidth',2, 'Color',[0 0 0])
plot(region_sizes, all_asphericity_PCA_mean(:,2), 'LineWidth',2, 'Color',[0.4 0.4 0.4])
plot(region_sizes, all_asphericity_PCA_mean(:,3), 'LineWidth',2, 'Color',[0.8 0.8 0.8])
plot(region_sizes, all_asphericity_PCA_mean(:,4), 'LineWidth',2, 'Color',[0.9290 0.6940 0.1250])
plot(region_sizes, all_asphericity_PCA_mean(:,5), 'LineWidth',2, 'Color', [0 0.4470 0.7410])
plot(region_sizes, all_asphericity_PCA_mean(:,6), 'LineWidth',2, 'Color', [0.4660 0.6740 0.1880])
xline((mean(asphericity_PCA_TADs_microadriaticum(:,1))+mean(asphericity_PCA_TADs_kawagutii(:,1)))/2,'--',{'Mean TAD Size'});
xlabel('Region Size [bp]', 'fontsize', 24)
ylabel('Asphericity (PCA)', 'fontsize', 24)
ax = gca;
ax.FontSize = 14;
legend('random regions CLC (16 discs)','random regions CLC (27 discs)','random regions CLC (54 discs)','random regions equilibrium globule', 'random regions S. microadriaticum', 'random regions S. kawagutii')
%legend('random regions CLC (16 discs)','random regions CLC (27 discs)','random regions CLC (54 discs)','random regions equilibrium globule', 'random regions S. microadriaticum', 'random regions S. kawagutii', 'whole chromosome S. microadriaticum', 'whole chromosome S. kawagutii', 'TADs S. microadriaticum', 'TADs S. kawagutii')
ylim([1 7.5])
legend boxoff
hold off




