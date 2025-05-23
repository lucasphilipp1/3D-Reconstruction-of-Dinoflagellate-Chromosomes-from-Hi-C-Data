% %Author: Lucas Philipp
% Fig S3 C,D,E,F in paper

% optomization of CSynth parameters using Hi-C data and multiplexed FISH data

% clc
% clear

% write human IMR90 chr 21 Hi-C data in CSynth compatible format
%%%
% chr21_Su_HiC = readtable('Hi-C_contacts_chromosome21.csv');
% chr21_Su_HiC=table2array(chr21_Su_HiC);
% %output for CSynth
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCSynth = zeros(nat_sum(size(chr21_Su_HiC,1)-1),3);
% count = 0;
% for i=2:1:size(chr21_Su_HiC,1)
%     for j=i:1:size(chr21_Su_HiC,1)
%         count = count + 1;
%         PCSynth(count,1)=chr21_Su_HiC(1,j);
%         PCSynth(count,2)=chr21_Su_HiC(i,1);
%         PCSynth(count,3)=chr21_Su_HiC(i,j);
%     end
% end
%
% writematrix(PCSynth,'Su_Cell_2020_chr21_CSynth_input.txt','Delimiter','tab')

%Comparison of CSynth Structures to FISH data for CSynth parameter optomization
%Using 651 FISH probe dataset: https://doi.org/10.1016/j.cell.2020.07.032

clc
clear

%determines how many FISH measurements to load, there are a lot!
rows=10000;

chromosome21_FISH_table = readtable("chromosome21.tsv", "FileType","text",'Delimiter', '\t');
total_measurements=size(chromosome21_FISH_table,1)
chromosome21_FISH = zeros(rows,6);

chromosome21_FISH(:,1) = chromosome21_FISH_table{1:rows,1}; %column 1: Z(nm)
chromosome21_FISH(:,2) =  chromosome21_FISH_table{1:rows,2}; %column 2: X(nm)
chromosome21_FISH(:,3) =  chromosome21_FISH_table{1:rows,3}; %column 3: Y(nm)
for i = 1:rows
    chromosome21_FISH(i,4) = str2double(extractBetween(chromosome21_FISH_table{i,4},':','-')); %column 4: base pair coordinate start
    chromosome21_FISH(i,5) = str2double(extractAfter(chromosome21_FISH_table{i,4},'-')); %column 5: base pair coordinate end
    i
end
chromosome21_FISH(:,6) =  chromosome21_FISH_table{1:rows,5}; %column 6: chromosome copy number

%remove NaNs
chromosome21_FISH=chromosome21_FISH(sum(isnan(chromosome21_FISH),2)==0,:);
%%%%%%%%%%%%

%chr21 how many probes?
num_probes = size(unique(chromosome21_FISH_table{1:rows,4}),1)
%chr21 how many cells?
cells = size(unique(chromosome21_FISH_table{1:rows,5}),1)

F = zeros(num_probes);
count = zeros(num_probes);

all_probes_start = unique(chromosome21_FISH(:,4));
all_probes_end = unique(chromosome21_FISH(:,5));

for i=1:1:cells
    %is there a measurement between probe j and probe k in cell i?
    indices = find(chromosome21_FISH(:,6) == i);
    probes = chromosome21_FISH(indices,4);

    for j=1:1:num_probes
        for k=1:1:num_probes
            if ismember(all_probes_start(j),probes) && ismember(all_probes_start(k),probes)
                pos1_rowindex = find(chromosome21_FISH(:,6) == i & chromosome21_FISH(:,4) == all_probes_start(j));
                pos2_rowindex = find(chromosome21_FISH(:,6) == i & chromosome21_FISH(:,4) == all_probes_start(k));

                count(j,k) = count(j,k) + 1;
                xdiff=abs(chromosome21_FISH(pos1_rowindex,2)-chromosome21_FISH(pos2_rowindex,2));
                ydiff=abs(chromosome21_FISH(pos1_rowindex,3)-chromosome21_FISH(pos2_rowindex,3));
                zdiff=abs(chromosome21_FISH(pos1_rowindex,1)-chromosome21_FISH(pos2_rowindex,1));

                F(j,k) = F(j,k) + sqrt(xdiff^2+ydiff^2+zdiff^2); %F, storage of FISH distances
            else
            end
        end
    end
    i
    
    if i == 1
    % read image and process it
    margins = [7 7 15 11]; % N S E W
    A = imread('/Users/lucasphilipp/Downloads/image.png');
    A = A(1+margins(1):end-margins(2),1+margins(4):end-margins(3),:);
    CT0 = permute(mean(im2double(A),1),[2 3 1]);
    CT0 = CT0([true; ~all(diff(CT0,1,1)==0,2)],:); % remove duplicate rows
    % CT0 is now a fixed-length color table
    % make it whatever length we want
    N = 32; % specify the number of colors in table
    na = size(CT0,1);
    CT = interp1(linspace(0,1,na),CT0,linspace(0,1,N));

    %Fig S3 C, visualization of FISH probes in a single cell
        chromosome21_FISH_plot=chromosome21_FISH(find(chromosome21_FISH(:,6)==1),:);
        figure
        scatter3(chromosome21_FISH_plot(:,1)-min(chromosome21_FISH_plot(:,1)),chromosome21_FISH_plot(:,2)-min(chromosome21_FISH_plot(:,2)),chromosome21_FISH_plot(:,3)-min(chromosome21_FISH_plot(:,3)),100,chromosome21_FISH_plot(:,4),'filled')
        colormap(CT)
        cb=colorbar;
        cb.Label.String = '[bp]';
        cb.FontSize = 18;
        caxis([0 46709983]); %hg38 assembly
        xlabel('x [nm]','FontSize', 24)
        ylabel('y [nm]','FontSize', 24)
        zlabel('z [nm]','FontSize', 24)
        view(45,15)
    end
end

%take cell average of FISH distances
F = F./count;
%diagonal of F should be 0

% figure
% fig=gcf;
% fig.Position(3:4)=[750,600];
% imagesc(F);
% colorbar
% xlabel('FISH probe #', 'fontsize', 24)
% ylabel('FISH probe #', 'fontsize', 24)
% title({'FISH probe separation'},{append('- ','chr21',' human IMR90 cells')}, 'FontSize', 24)
% cb=colorbar;
% cb.Label.String = 'Cell-averaged probe separation [nm]';
% cb.FontSize = 18;

%load CSynth structures generated with different parameters
CSynth_structures={
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_20_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_40_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_60_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_80_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_100_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_20_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_40_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_60_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_80_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_100_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_20_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_40_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_60_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_80_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_100_PP_-1.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_20_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_40_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_60_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_80_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_100_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_20_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_40_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_60_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_80_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_100_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_20_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_40_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_60_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_80_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_100_PP_-2.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_20_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_40_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_60_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_80_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_100_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_20_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_40_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_60_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_80_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_100_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_20_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_40_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_60_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_80_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_100_PP_-3.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_20_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_40_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_60_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_80_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_0_CF_100_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_20_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_40_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_60_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_80_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-1_CF_100_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_20_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_40_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_60_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_80_PP_-4.xyz' ...
    'Su_Cell_2020_chr21_CSynth_output_SP_-2_CF_100_PP_-4.xyz' ...
    };

avg_RE=zeros(size(CSynth_structures,2),1); %average relative error
CSynth_chr_21 = tdfread(CSynth_structures{1});
names_chr21 = fieldnames(CSynth_chr_21);
ALL_CSynth_structures=zeros(size(getfield(CSynth_chr_21,names_chr21{1}),1),4,size(CSynth_structures,2));

%CSynth predicted 3D structure
%column 1: bp, column 2,3,4: x y z position
for q = 1:1:size(CSynth_structures,2)
    CSynth_chr_21 = tdfread(CSynth_structures{q});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    names_chr21 = fieldnames(CSynth_chr_21);
    HiC = zeros(size(getfield(CSynth_chr_21,names_chr21{1}),1),4);
    HiC(:,1) = getfield(CSynth_chr_21,names_chr21{1}); %column 1: bp
    HiC(:,2) = getfield(CSynth_chr_21,names_chr21{2}); %column 2: x position, arb units
    HiC(:,3) = getfield(CSynth_chr_21,names_chr21{3}); %column 3: y position, arb units
    HiC(:,4) = getfield(CSynth_chr_21,names_chr21{4}); %column 4: z position, arb units

    H = zeros(num_probes); %H, storage of CSynth-predicted distances from Hi-C

    %FISH probes are not genomically uniformly distributed
    %map HiC coordinates to individual probes
    for i=1:1:num_probes
        for j=1:1:num_probes

            keep_i = find(HiC(:,1) >= all_probes_start(i) & HiC(:,1) <= all_probes_end(i));
            keep_j = find(HiC(:,1) >= all_probes_start(j) & HiC(:,1) <= all_probes_end(j));

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

    %rescale CSynth coordinates isotropically to best-match with microscopy
    HiC(:,2) = HiC(:,2)*scale_Factor;
    HiC(:,3) = HiC(:,3)*scale_Factor;
    HiC(:,4) = HiC(:,4)*scale_Factor;

    %save rescaled structure
    ALL_CSynth_structures(:,:,q)=HiC;

    %recompute CSynth probe distances will rescaled axes
    H = zeros(num_probes);
    for i=1:1:num_probes
        for j=1:1:num_probes

            keep_i = find(HiC(:,1) >= all_probes_start(i) & HiC(:,1) <= all_probes_end(i));
            keep_j = find(HiC(:,1) >= all_probes_start(j) & HiC(:,1) <= all_probes_end(j));

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

    % figure
    % fig=gcf;
    % fig.Position(3:4)=[750,600];
    % imagesc(H);
    % colorbar
    % xlabel('FISH probe #', 'fontsize', 24)
    % ylabel('FISH probe #', 'fontsize', 24)
    % title({'CSynth predicted probe separation using HiC data'}, {append('- ','chr21',' human IMR90 cells')}, 'FontSize', 24)
    % cb=colorbar;
    % cb.Label.String = 'Probe separation [CSynth units]';
    % cb.FontSize = 18;

    %relative error, RE
    %RE=|H_ij - F_ij|/F_ij;
    %F_ij FISH distance, H_ij HiC distance, ij genomic loci
    RE = abs(H-F)./F;

    RE(idiag(size(RE),-1)) = NaN(size(RE,1)-1,1);
    RE(idiag(size(RE),1)) = NaN(size(RE,1)-1,1);

    RE_vec = RE(:);
    F_vec = F(:);
    H_vec = H(:);

    idx=find(F_vec == 0);
    F_vec(idx) = [];
    H_vec(idx) = [];

    %NaNs due to diagonal ommited from average
    RE_vec=RE_vec(~isinf(RE_vec));
    RE_vec=RE_vec(~isnan(RE_vec));

    %we want to minimize the average relative error as a function of CSynth parameters
    %disp('Average relative error')
    avg_RE(q)=median(RE_vec)
    q

    %diagonal set to zero for plotting purposes
    RE(isinf(RE))=0;

    %Fig S3 D & F
    % figure
    % fig=gcf;
    % fig.Position(3:4)=[750,600];
    % imagesc(RE);
    % xlabel('FISH probe #', 'fontsize', 24)
    % ylabel('FISH probe #', 'fontsize', 24)
    % title({'Relative Error'},{append('- ',CSynth_structures{q},' human IMR90 cells')}, 'FontSize', 24, 'Interpreter', 'none')
    % cb=colorbar;
    % cb.Label.String = 'Relative Error';
    % cb.FontSize = 18;
    % caxis([0 1]);

    %are larger distortions happening at larger probe separations?
    % figure
    % hist3([F_vec, RE_vec],[100 100],'CdataMode','auto')
    % view(2)
    % xlabel('FISH probe separation [nm]', 'fontsize', 24)
    % ylabel('Relative Error = |H_{ij} - F_{ij}|/F_{ij}', 'fontsize', 24)
    % title({
    %     ['Are larger distortions happening']
    %     ['at larger probe separations?']
    %     }, 'fontsize', 24);
    % subtitle(append('- ',CSynth_structures{q},' human IMR90 cells'), 'Interpreter', 'none')
    % xlim([min(F_vec) max(F_vec)])
    % ylim([0 max(RE_vec)])
end

% plot all CSynth structures simultaneously, Fig S3 E
%%
figure
hold on
for q = 1:1:size(CSynth_structures,2)
    plot3(ALL_CSynth_structures(:,2,q),ALL_CSynth_structures(:,3,q),ALL_CSynth_structures(:,4,q),'-')
end
xlabel('x', 'fontsize', 18)
ylabel('y', 'fontsize', 18)
zlabel('z', 'fontsize', 18)
title('CSynth 3D Prediction - IMR90 HiC chr 21')
grid on
temp_x = ALL_CSynth_structures(:,2,:);
temp_y = ALL_CSynth_structures(:,3,:);
temp_z = ALL_CSynth_structures(:,4,:);
xlim([min(temp_x(:)) max(temp_x(:))])
ylim([min(temp_y(:)) max(temp_y(:))])
zlim([min(temp_z(:)) max(temp_z(:))])
view(-30,15)
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])
axis equal

param_avg_RE=zeros(size(CSynth_structures,2),4);
param_avg_RE(:,1) = avg_RE;
param_avg_RE(:,2) = [0 0 0 0 0 -1 -1 -1 -1 -1 -2 -2 -2 -2 -2 0 0 0 0 0 -1 -1 -1 -1 -1 -2 -2 -2 -2 -2 0 0 0 0 0 -1 -1 -1 -1 -1 -2 -2 -2 -2 -2 0 0 0 0 0 -1 -1 -1 -1 -1 -2 -2 -2 -2 -2]; %SP, spring power
param_avg_RE(:,3) = [20 40 60 80 100 20 40 60 80 100 20 40 60 80 100 20 40 60 80 100 20 40 60 80 100 20 40 60 80 100 20 40 60 80 100 20 40 60 80 100 20 40 60 80 100 20 40 60 80 100 20 40 60 80 100 20 40 60 80 100]; %CF, contact force
param_avg_RE(:,4) = [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4]; %PP, pushapart power

%visualization of relative error v.s. CSynth parameters
figure
scatter3(param_avg_RE(:,2),param_avg_RE(:,3),param_avg_RE(:,4),500,param_avg_RE(:,1),'filled')
ax = gca;
view(-30,15)
xlabel('spring power','fontsize', 18)
ylabel('contact force', 'fontsize', 18)
zlabel('pushapart power', 'fontsize', 18)
cb = colorbar;
cb.Label.String = 'Relative Error';
cb.FontSize = 14;
%caxis([0 1]);

function sum = nat_sum(x)
sum=0;
for i=1:x
    sum=sum+i;
end
end

function [I J] = idiag(sz, k)
if isscalar(sz)
    sz = [sz sz];
end
m=sz(1);
n=sz(2);
if nargin<2
    k=0;
end
I = (1-min(k,0):min(m,n-k)).';
J = I+k;
if nargout<2
    I = sub2ind([m n], I, J);
end
end

