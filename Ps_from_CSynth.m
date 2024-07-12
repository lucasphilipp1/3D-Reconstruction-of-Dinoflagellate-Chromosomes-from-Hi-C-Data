clc
clear

num_chroms = 75; %omit small HiC scaffolds that probably aren't complete chromosomes

resolution = 5000; %HiC resolution in bp
s_avg = cell(1, num_chroms); %cell array for separation
Ps_avg = cell(1, num_chroms); %cell array for contact probability
cont_count_cell = cell(1, num_chroms);

adj_monomer_sep_list=[];
chr_sizes=[];

for i = 1:num_chroms
    chromosome = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i.xyz',i));
    %chromosome = importdata(sprintf('symbiodinium_microadriaticum_chr%i_3D.xyz',i));
    for j = 1:1:size(chromosome,1)-1
        adj_monomer_sep_list = [adj_monomer_sep_list; pdist(chromosome(j:j+1,2:4))];
    end
    chr_sizes = [chr_sizes; size(chromosome,1)];
end
distance_threshold = mean(adj_monomer_sep_list)*3; %choose the smallest distance threshold where no shot noise is observed
min_chr_size = min(chr_sizes);
max_chr_size = max(chr_sizes);

for i = 1:num_chroms
    chromosome = importdata(sprintf('s_kawagutii_V3_HiC_scaffold_%i.xyz',i));
    %chromosome = importdata(sprintf('symbiodinium_microadriaticum_chr%i_3D.xyz',i));
    
    D = pdist(chromosome(:,2:4));
    D = squareform(D);
    cont_count_mat = D;
    cont_count_mat(cont_count_mat<distance_threshold) = 0;
    cont_count_mat(cont_count_mat>distance_threshold) = 1;
    cont_count_mat(:) = ~cont_count_mat;

    cont_count=zeros(1,size(chromosome(:,2:4),1)-1);
    for j=1:1:size(chromosome(:,2:4),1)-1
        cont_count(j) = sum(diag(cont_count_mat,j));
    end
    cont_count_cell{i}=cont_count';
end

bins=[4];
while bins(end) < max_chr_size-100
    bins=[bins; bins(end)*1.05];
end
bins = unique(round(bins));
bins(end)=[];

Ps = size(bins,1)-1;

for i = 1:size(bins,1)-1
    s=bins(i);
    while s<bins(i+1)
        s = [s; s(end)+1];
    end
    s(end)=[];

    contacts = 0;
    possible_mon_pairs=0;
    for j = 1:num_chroms
        if max(s)<size(cont_count_cell{1,j},1)
            contacts = contacts+sum(cont_count_cell{1,j}(s));
            for k = 1:size(s)
                possible_mon_pairs = possible_mon_pairs + chr_sizes(j)-s(k);
            end
        end
    end
    Ps(i)=contacts./possible_mon_pairs;
end

figure
hold on
ax=gca;
A=4; %vertical offset for dotted line
xA=linspace(bins(1)*resolution,8*10^6,100);
plot(xA,A*xA.^(-0.5), '--') %mitotic chromosome exponent
plot(bins(1:end-15).*resolution,Ps(1:end-14),'LineWidth', 2)
%legendStrings = "d_{cutoff} = " + string(distance_threshold) + "um";
%legend(legendStrings, 'Location', 'southwest')
xlabel('s [bp]', 'fontsize', 24)
ylabel('P(s)', 'fontsize', 24)
set(ax,'xScale', 'log')
set(ax,'YScale', 'log')
ax = gca;
ax.FontSize = 24; 
xlim([10^4 2*10^7])
ylim([10^-6 10^0])
grid on