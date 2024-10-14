%Fig 3D & S9 in paper
%visualization of spatial variation in transcription levels using a
%cylindrical coorindate system

clc
clear

num_chroms = 75;
bins = 20;

%symbiodinium microadriaticum
all_expression_agg = [];

telomeres = [];

for i=1:num_chroms

    chromosome = importdata(sprintf('symbiodinium_microadriaticum_chr%i_3D.xyz',i));
    RNAseq = importdata(sprintf('s_microadriaticum_chr%i_pilon_strandsum_TPM.bed',i));

    chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
    chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));
    chromosome(:,4) = chromosome(:,4) - mean(chromosome(:,4));

    HiC_resolution = chromosome(2,1)-chromosome(1,1);

    genes_start_end_TPM = [RNAseq.data(:,1) RNAseq.data(:,2) RNAseq.data(:,8)];

    [coeff,score] = pca(chromosome(:,2:4));
    chromosome_PCA = score;

    %get xyz coordinates for every base pair by linearly interpolating between the model
    chromosome_PCA_interpolated = zeros(chromosome(end,1),4);
    chromosome_PCA_interpolated(:,1) = linspace(1,chromosome(end,1),chromosome(end,1))';
    for j = 1:1:size(chromosome_PCA,1)-1
        chromosome_PCA_interpolated(HiC_resolution*(j-1)+1:HiC_resolution*j,2:4) = [linspace(chromosome_PCA(j,1),chromosome_PCA(j+1,1),HiC_resolution)' linspace(chromosome_PCA(j,2),chromosome_PCA(j+1,2),HiC_resolution)' linspace(chromosome_PCA(j,3),chromosome_PCA(j+1,3),HiC_resolution)'];
    end

    %some sequence at end of chromosome is truncated because monomers are 1/5kbp
    while genes_start_end_TPM(end,2) > chromosome_PCA_interpolated(end,1)
        genes_start_end_TPM(end,:)=[];
    end

    genes_center_xyz = [];

    %store coordinates for each base pair in each gene
    for k = 1:1:size(genes_start_end_TPM,1)
        genes_center_xyz = [genes_center_xyz; mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),2)) mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),3)) mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),4))];
    end

    %column 1 is PC1, the chromosome long axis
    %cart2pol does x,y,z -> theta,rho,z. z is unchanged
    [theta_TPM, rho_TPM, z_TPM] = cart2pol(genes_center_xyz(:,2),genes_center_xyz(:,3),genes_center_xyz(:,1));

    %column 1 is PC1, the chromosome long axis
    [theta_chromosome, rho_chromosome, z_chromosome] = cart2pol(chromosome_PCA(:,2),chromosome_PCA(:,3),chromosome_PCA(:,1));

    temp_telomeres = [chromosome_PCA(1,:); chromosome_PCA(end,:)];

    %rescale chromosome dimensions
    all_expression_agg = [all_expression_agg; rho_TPM./max(rho_chromosome) z_TPM./max(abs(z_chromosome)) abs(genes_start_end_TPM(:,3)) ones(size(genes_start_end_TPM,1),1).*i];
    %telomeres were defined as = 1% from the start/end of the chromosome which is on average ~ 100kbp.
    telomeres = [telomeres; rho_chromosome(round(size(rho_chromosome,1)*0.01),:)./max(rho_chromosome) z_chromosome(round(size(rho_chromosome,1)*0.01),:)./max(abs(z_chromosome))];
    telomeres = [telomeres; rho_chromosome(end-round(size(rho_chromosome,1)*0.01),:)./max(rho_chromosome) z_chromosome(end-round(size(rho_chromosome,1)*0.01),:)./max(abs(z_chromosome))];
    i
end

%without telomeres
figure
[N,c] = hist3(all_expression_agg(:,1:2),[bins bins]);
N=flipud(N);
c{1}=fliplr(c{1});
for i = 1:1:size(N,1)
    N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
end
buffer_length=7;
Nbuffer = cat(1,zeros(buffer_length,20),N); %add a buffer of zeros to include more of chromosome periphery in plot
imagesc(Nbuffer)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
% hold on
% plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1+buffer_length,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k', 'MarkerSize', 4)
% hold off
colormap("sky")
cb = colorbar;
cb.Label.String = 'Normalized Gene Density';
cb.FontSize = 24;
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('cylidrical axis','FontSize', 24)
ylabel('radial axis','FontSize', 24)
t=get(cb,'Limits');
set(cb,'Ticks',[t(1),t(2)])
set(cb,'XTickLabel',{'Low','High',});
%xlim([2.5 20.5])
pbaspect([1.3 1 1])

%with telomeres
figure
[N,c] = hist3(all_expression_agg(:,1:2),[bins bins]);
N=flipud(N);
c{1}=fliplr(c{1});
for i = 1:1:size(N,1)
    N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
end
buffer_length=7;
Nbuffer = cat(1,zeros(buffer_length,20),N); %add a buffer of zeros to include more of chromosome periphery in plot
imagesc(Nbuffer)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
hold on
plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1+buffer_length,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k', 'MarkerSize', 4)
hold off
colormap("sky")
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
pbaspect([1.3 1 1])

%symbiodinium kawagutii
i=0;
all_expression_agg = [];

telomeres = [];

for i=1:num_chroms

    chromosome = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i.xyz',i));
    RNAseq = importdata(sprintf('s_kawagutii_chr%i_strandsum_TPM.bed',i));

    chromosome(:,2) = chromosome(:,2) - mean(chromosome(:,2));
    chromosome(:,3) = chromosome(:,3) - mean(chromosome(:,3));
    chromosome(:,4) = chromosome(:,4) - mean(chromosome(:,4));

    HiC_resolution = chromosome(2,1)-chromosome(1,1);

    genes_start_end_TPM = [RNAseq.data(:,1) RNAseq.data(:,2) RNAseq.data(:,8)];

    [coeff,score] = pca(chromosome(:,2:4));
    chromosome_PCA = score;

    %get xyz coordinates for every base pair by linearly interpolating between the model
    chromosome_PCA_interpolated = zeros(chromosome(end,1),4);
    chromosome_PCA_interpolated(:,1) = linspace(1,chromosome(end,1),chromosome(end,1))';
    for j = 1:1:size(chromosome_PCA,1)-1
        chromosome_PCA_interpolated(HiC_resolution*(j-1)+1:HiC_resolution*j,2:4) = [linspace(chromosome_PCA(j,1),chromosome_PCA(j+1,1),HiC_resolution)' linspace(chromosome_PCA(j,2),chromosome_PCA(j+1,2),HiC_resolution)' linspace(chromosome_PCA(j,3),chromosome_PCA(j+1,3),HiC_resolution)'];
    end

    while genes_start_end_TPM(end,2) > chromosome_PCA_interpolated(end,1)
        genes_start_end_TPM(end,:)=[];
    end

    genes_center_xyz = [];

    %store coordinates for each base pair in each gene
    for k = 1:1:size(genes_start_end_TPM,1)
        genes_center_xyz = [genes_center_xyz; mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),2)) mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),3)) mean(chromosome_PCA_interpolated(genes_start_end_TPM(k,1):genes_start_end_TPM(k,2),4))];
    end

    %column 1 is PC1, the chromosome long axis
    %cart2pol does x,y,z -> theta,rho,z. z is unchanged
    [theta_TPM, rho_TPM, z_TPM] = cart2pol(genes_center_xyz(:,2),genes_center_xyz(:,3),genes_center_xyz(:,1));

    %column 1 is PC1, the chromosome long axis
    [theta_chromosome, rho_chromosome, z_chromosome] = cart2pol(chromosome_PCA(:,2),chromosome_PCA(:,3),chromosome_PCA(:,1));

    temp_telomeres = [chromosome_PCA(1,:); chromosome_PCA(end,:)];

    %rescale chromosome dimensions
    all_expression_agg = [all_expression_agg; rho_TPM./max(rho_chromosome) z_TPM./max(abs(z_chromosome)) abs(genes_start_end_TPM(:,3)) ones(size(genes_start_end_TPM,1),1).*i];
    %telomeres were defined as = 1% from the start/end of the chromosome which is on average ~ 100kbp.
    telomeres = [telomeres; rho_chromosome(round(size(rho_chromosome,1)*0.01),:)./max(rho_chromosome) z_chromosome(round(size(rho_chromosome,1)*0.01),:)./max(abs(z_chromosome))];
    telomeres = [telomeres; rho_chromosome(end-round(size(rho_chromosome,1)*0.01),:)./max(rho_chromosome) z_chromosome(end-round(size(rho_chromosome,1)*0.01),:)./max(abs(z_chromosome))];
    i
end

%without telomeres
figure
[N,c] = hist3(all_expression_agg(:,1:2),[bins bins]);
N=flipud(N);
c{1}=fliplr(c{1});
for i = 1:1:size(N,1)
    N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
end
buffer_length=7;
Nbuffer = cat(1,zeros(buffer_length,20),N); %add a buffer of zeros to include more of chromosome periphery in plot
imagesc(Nbuffer)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
% hold on
% plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1+buffer_length,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k', 'MarkerSize', 4)
% hold off
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
pbaspect([1.3 1 1])

%with telomeres
figure
[N,c] = hist3(all_expression_agg(:,1:2),[bins bins]);
N=flipud(N);
c{1}=fliplr(c{1});
for i = 1:1:size(N,1)
    N(:,i)=N(:,i)./c{1}'; %cylidrical jacobian correction
end
buffer_length=7;
Nbuffer = cat(1,zeros(buffer_length,20),N); %add a buffer of zeros to include more of chromosome periphery in plot
imagesc(Nbuffer)
%r coordinate mapping: [0 1]->[bins 1]
%z coordinate mapping: [-1 1]->[1 bins]
hold on
plot((telomeres(:,2).*((bins-1)/2))+(bins-1)/2+1,((-1.*telomeres(:,1))+1).*(bins-1)+1+buffer_length,'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k', 'MarkerSize', 4)
hold off
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
pbaspect([1.3 1 1])

%for green colormap
%https://github.com/JonathanRob/GeneSetAnalysisMatlab/blob/master/custom_cmap.m
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


