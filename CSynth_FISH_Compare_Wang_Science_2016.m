clc
clear

%Comparison of CSynth Structures to FISH data

CSynth_chr_20 = tdfread('IMR90_5kbp_chr20_MAPQG0_RAW_3D.xyz');
names_chr20 = fieldnames(CSynth_chr_20);
chr20_HiC = zeros(size(getfield(CSynth_chr_20,names_chr20{1}),1),4); %column 1: bp, column 2,3,4: x y z position
chr20_HiC(:,1) = getfield(CSynth_chr_20,names_chr20{1});
chr20_HiC(:,2) = getfield(CSynth_chr_20,names_chr20{2});
chr20_HiC(:,3) = getfield(CSynth_chr_20,names_chr20{3});
chr20_HiC(:,4) = getfield(CSynth_chr_20,names_chr20{4});

figure
plot3(chr20_HiC(:,2),chr20_HiC(:,3),chr20_HiC(:,4),'-')
xlabel('x')
ylabel('y')
zlabel('z')
title('CSynth 3D Prediction - IMR90 HiC chr 20')
grid on
axis equal

CSynth_chr_21 = tdfread('IMR90_5kbp_chr21_MAPQG0_RAW_3D.xyz');
names_chr21 = fieldnames(CSynth_chr_21);
chr21_HiC = zeros(size(getfield(CSynth_chr_21,names_chr21{1}),1),4); %column 1: bp, column 2,3,4: x y z position
chr21_HiC(:,1) = getfield(CSynth_chr_21,names_chr21{1});
chr21_HiC(:,2) = getfield(CSynth_chr_21,names_chr21{2});
chr21_HiC(:,3) = getfield(CSynth_chr_21,names_chr21{3});
chr21_HiC(:,4) = getfield(CSynth_chr_21,names_chr21{4});

figure
plot3(chr21_HiC(:,2),chr21_HiC(:,3),chr21_HiC(:,4),'-')
xlabel('x')
ylabel('y')
zlabel('z')
title('CSynth 3D Prediction - IMR90 HiC chr 21')
grid on
axis equal

CSynth_chr_22 = tdfread('IMR90_5kbp_chr22_MAPQG0_RAW_3D.xyz');
names_chr22 = fieldnames(CSynth_chr_22);
chr22_HiC = zeros(size(getfield(CSynth_chr_22,names_chr22{1}),1),4); %column 1: bp, column 2,3,4: x y z position
chr22_HiC(:,1) = getfield(CSynth_chr_22,names_chr22{1});
chr22_HiC(:,2) = getfield(CSynth_chr_22,names_chr22{2});
chr22_HiC(:,3) = getfield(CSynth_chr_22,names_chr22{3});
chr22_HiC(:,4) = getfield(CSynth_chr_22,names_chr22{4});

figure
plot3(chr22_HiC(:,2),chr22_HiC(:,3),chr22_HiC(:,4),'-')
xlabel('x')
ylabel('y')
zlabel('z')
title('CSynth 3D Prediction - IMR90 HiC chr 22')
grid on
axis equal

%bp file
%ID# of FISH probe, Start gemomic coordinate [bp], End genomic coordinate [bp]

%xyz file
%Cell #, ID# of FISH probe, x [um], y [um], z [um]
%must average data across multiple chromosomes

%chr20bp = readmatrix('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/Positive_Control_ Wang_ Science_2016/IMR90/FISH/chr20bp.csv');
%chr20xyz = readmatrix('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/Positive_Control_ Wang_ Science_2016/IMR90/FISH/chr20xyz.csv');
chr21bp = readmatrix('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/Positive_Control_ Wang_ Science_2016/IMR90/FISH/chr21bp.csv');
chr21xyz = readmatrix('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/Positive_Control_ Wang_ Science_2016/IMR90/FISH/chr21xyz.csv');
%chr22bp = readmatrix('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/Positive_Control_ Wang_ Science_2016/IMR90/FISH/chr22bp.csv');
%chr22xyz = readmatrix('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/Positive_Control_ Wang_ Science_2016/IMR90/FISH/chr22xyz.csv');

%remove NaNs
%chr20xyz=chr20xyz(sum(isnan(chr20xyz),2)==0,:);
chr21xyz=chr21xyz(sum(isnan(chr21xyz),2)==0,:);
%chr22xyz=chr22xyz(sum(isnan(chr22xyz),2)==0,:);

%%%%%%%%%%%%

for_figure_caption={'chr20','chr21','chr22'};

%for q=1:3
q=2;
    if q==1
        %chr20 how many probes?
        num_probes = size(chr20bp,1);
        %chr20 how many cells?
        cells = size(unique(chr20xyz(:,1)),1);

        bp=chr20bp;
        xyz=chr20xyz;
        HiC=chr20_HiC;
    elseif q==2
        %chr21 how many probes?
        num_probes = size(chr21bp,1);
        %chr21 how many cells?
        cells = size(unique(chr21xyz(:,1)),1);

        bp=chr21bp;
        xyz=chr21xyz;
        HiC=chr21_HiC;
    else
        %chr21 how many probes?
        num_probes = size(chr22bp,1);
        %chr21 how many cells?
        cells = size(unique(chr22xyz(:,1)),1);

        bp=chr22bp;
        xyz=chr22xyz;
        HiC=chr22_HiC;
    end

    F = zeros(num_probes);
    count = zeros(num_probes);

    for i=1:1:cells
        %is there a measurement between probe j and probe k in cell i?
        indices = find(xyz(:,1) == i);
        probes = xyz(indices,2);
        for j=1:1:num_probes
            for k=1:1:num_probes
                if ismember(j,probes) && ismember(k,probes)
                    pos1_rowindex = find(xyz(:,1) == i & xyz(:,2) == j);
                    pos2_rowindex = find(xyz(:,1) == i & xyz(:,2) == k);

                    count(j,k) = count(j,k) + 1;
                    xdiff=abs(xyz(pos1_rowindex,3)-xyz(pos2_rowindex,3));
                    ydiff=abs(xyz(pos1_rowindex,4)-xyz(pos2_rowindex,4));
                    zdiff=abs(xyz(pos1_rowindex,5)-xyz(pos2_rowindex,5));

                    F(j,k) = F(j,k) + sqrt(xdiff^2+ydiff^2+zdiff^2);
                else
                end
            end
        end
    end

    %take cell average of FISH distances
    F = F./count;
    %diagonal of F should be 0

    H = zeros(num_probes);

    for i=1:1:num_probes
        for j=1:1:num_probes

            keep_i = find(HiC(:,1) >= bp(i,2) & HiC(:,1) <= bp(i,3));
            keep_j = find(HiC(:,1) >= bp(j,2) & HiC(:,1) <= bp(j,3));

            for k = 1:1:size(keep_i,1)
                for l = 1:1:size(keep_j,1)
                    xdiff=abs(HiC(keep_i(k),2)-HiC(keep_j(l),2));
                    ydiff=abs(HiC(keep_i(k),3)-HiC(keep_j(l),3));
                    zdiff=abs(HiC(keep_i(k),4)-HiC(keep_j(l),4));
                    H(i,j) = H(i,j) + sqrt(xdiff^2+ydiff^2+zdiff^2);
                end
            end

            H(i,j) = H(i,j)./(size(keep_i,1)*size(keep_j,1)); %divides by zero if HiC data is lower resolution than FISH data

        end
    end

    figure
    fig=gcf;
    fig.Position(3:4)=[750,600];
    imagesc(H);
    colorbar
    xlabel('FISH probe #', 'fontsize', 24)
    ylabel('FISH probe #', 'fontsize', 24)
    title({'CSynth predicted probe separation using HiC data'}, {append('- ',for_figure_caption{q},' human IMR90 cells')}, 'FontSize', 24)
    cb=colorbar;
    cb.Label.String = 'Probe separation [CSynth units]';
    cb.FontSize = 18;

    F_vec = F(:);
    H_vec = H(:);
    idx=find(F_vec == 0);
    F_vec(idx) = [];
    H_vec(idx) = [];

    %mean_separate performs best
    mean_together=mean(F_vec./H_vec);
    mean_seperate=mean(F_vec)/mean(H_vec);
    median_together=median(F_vec./H_vec);
    median_seperate=median(F_vec)/median(H_vec);

    scale_Factor = mean_seperate;

    %rescale CSynth coordinates, all axes equally
    HiC(:,2) = HiC(:,2)*scale_Factor;
    HiC(:,3) = HiC(:,3)*scale_Factor;
    HiC(:,4) = HiC(:,4)*scale_Factor;

    %recompute CSynth probe distances will rescaled axes
    H = zeros(num_probes);

    for i=1:1:num_probes
        for j=1:1:num_probes

            keep_i = find(HiC(:,1) >= bp(i,2) & HiC(:,1) <= bp(i,3));
            keep_j = find(HiC(:,1) >= bp(j,2) & HiC(:,1) <= bp(j,3));

            for k = 1:1:size(keep_i,1)
                for l = 1:1:size(keep_j,1)
                    xdiff=abs(HiC(keep_i(k),2)-HiC(keep_j(l),2));
                    ydiff=abs(HiC(keep_i(k),3)-HiC(keep_j(l),3));
                    zdiff=abs(HiC(keep_i(k),4)-HiC(keep_j(l),4));
                    H(i,j) = H(i,j) + sqrt(xdiff^2+ydiff^2+zdiff^2);
                end
            end

            H(i,j) = H(i,j)./(size(keep_i,1)*size(keep_j,1)); %divides by zero if HiC data is lower resolution than FISH data

        end
    end

    figure
    fig=gcf;
    fig.Position(3:4)=[750,600];
    imagesc(F);
    colorbar
    xlabel('FISH probe #', 'fontsize', 24)
    ylabel('FISH probe #', 'fontsize', 24)
    title({'FISH probe separation'},{append('- ',for_figure_caption{q},' human IMR90 cells')}, 'FontSize', 24)
    cb=colorbar;
    cb.Label.String = 'Cell-averaged probe separation [um]';
    cb.FontSize = 18;

    %relative error, RE
    %RE=|H_ij - F_ij|/F_ij;
    %F_ij FISH distance, H_ij HiC distance, ij genomic loci
    RE = abs(H-F)./F;

    RE_vec = RE(:);
    F_vec = F(:);
    H_vec = H(:);

    idx=find(F_vec == 0);
    F_vec(idx) = [];
    H_vec(idx) = [];

    %NaNs due to diagonal ommited from average
    RE_vec=RE_vec(~isinf(RE_vec));

    %we want to minimize the average relative error as a function of CSynth parameters
    mean(RE_vec)

    %diagonal set to zero for plotting purposes
    RE(isinf(RE))=0;

    figure
    fig=gcf;
    fig.Position(3:4)=[750,600];
    imagesc(RE);
    xlabel('FISH probe #', 'fontsize', 24)
    ylabel('FISH probe #', 'fontsize', 24)
    title({'Relative Error'},{append('- ',for_figure_caption{q},' human IMR90 cells')}, 'FontSize', 24)
    cb=colorbar;
    cb.Label.String = 'abs(H_{ij}-F_{ij})/F_{ij}';
    cb.FontSize = 18;

    %are larger distortions happening at larger probe separations?
    figure
    scatter(F_vec, RE_vec,'filled')
    xlabel('FISH probe separation [um]', 'fontsize', 24)
    ylabel('Relative Error = |H_{ij} - F_{ij}|/F_{ij}', 'fontsize', 24)
    title({
        ['Are larger distortions happening']
        ['at larger probe separations?']
        }, 'fontsize', 24);
    subtitle(append('- ',for_figure_caption{q},' human IMR90 cells'))

    %pearson correlation between H and F matrix, ignoring diagonal entries
    %this does not depend on scaleFactor
    pearson=corrcoef(H_vec,F_vec);
    pearson=pearson(1,2)

    %this artificially increases the correlation by forcing diagonal
    %entries to be the same
    % H = H - diag(diag(H));
    % F = F - diag(diag(F));
    % 
    % matcorr=corrcoef(H,F)
%end
