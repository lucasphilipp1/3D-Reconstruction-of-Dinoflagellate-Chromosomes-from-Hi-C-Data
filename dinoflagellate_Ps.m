clc
clear

%f = fullfile('myfolder','mysubfolder','myfile.m')

resolution = 1000;

microadriaticum_coccoid_Ps = {};
microadriaticum_coccoid_s = {};

source_dir = ('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/hic contact probabilities microadriaticum/coccoid/');
N_microadriaticum_coccoid=size(dir([source_dir, '/*.txt']),1);
file_suffix = cell(N_microadriaticum_coccoid,1);

N_microadriaticum_coccoid=94;

for i=1:1:N_microadriaticum_coccoid
    file_suffix{i} = append(append('symbiodinium_microadriaticum_coccoid_chr',num2str(i)),'.txt');
    microadriaticum_coccoid_struct=tdfread(append('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/hic contact probabilities microadriaticum/coccoid/',file_suffix{i}));
    microadriaticum_coccoid_names = fieldnames(microadriaticum_coccoid_struct);

    col1 = getfield(microadriaticum_coccoid_struct,microadriaticum_coccoid_names{1});
    col2 = getfield(microadriaticum_coccoid_struct,microadriaticum_coccoid_names{2});
    col3 = getfield(microadriaticum_coccoid_struct,microadriaticum_coccoid_names{3});

    col_all = [col1 col2 col3];

    P = zeros(max(col1)/resolution);
    for j = 1:1:size(col_all,1)
        P(col1(j)/resolution,col2(j)/resolution) = col3(j);
        P(col2(j)/resolution,col1(j)/resolution) = col3(j);
        j
    end

    [s Ps]=contact_probability_HiC(P,resolution);
    
    s = s(2:end); 
    Ps = Ps(2:end);
    Ps=Ps/Ps(1); %P(min(s))=1
    
    writematrix([s round(Ps,3,"significant")],append(append('Ps_symbiodinium_microadriaticum_coccoid_chr',num2str(i)),'.txt'),'Delimiter','tab')

    figure
    hold on
    ax=gca;
    set(ax,'xScale', 'log')
    set(ax,'YScale', 'log')
    plot(s,Ps)
    xlabel('s')
    ylabel('P(s)')
    xlim([10^(3) 10^(8)])
    ylim([10^(-6) 10^(-1)])

    saveas(gcf,append('Ps_symbiodinium_microadriaticum_coccoid_chr',num2str(i)),'png')

    i
end

resolution = 5000;

bminutum_Ps = {};
bminutum_s = {};

source_dir = ('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/hic contact probabilities breviolum/');
N_bminutum=size(dir([source_dir, '/*.txt']),1);
file_suffix = cell(N_bminutum,1);

N_bminutum=100;

for i=1:1:N_bminutum
    file_suffix{i} = append(append('bminutum_pseudochromosome_',num2str(i)),'.txt');
    bminutum_struct=tdfread(append('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/hic contact probabilities breviolum/',file_suffix{i}));
    bminutum_names = fieldnames(bminutum_struct);

    col1 = getfield(bminutum_struct,bminutum_names{1});
    col2 = getfield(bminutum_struct,bminutum_names{2});
    col3 = getfield(bminutum_struct,bminutum_names{3});

    col_all = [col1 col2 col3];

    P = zeros(max(col1)/resolution);
    for j = 1:1:size(col_all,1)
        P(col1(j)/resolution+1,col2(j)/resolution+1) = col3(j);
        P(col2(j)/resolution+1,col1(j)/resolution+1) = col3(j);
        j
    end

    [s Ps]=contact_probability_HiC(P,resolution);

    s = s(2:end); 
    Ps = Ps(2:end);
    Ps=Ps/Ps(1); %P(min(s))=1

    writematrix([s Ps],append(append('Ps_bminutum_chr',num2str(i)),'.txt'),'Delimiter','tab')

    figure
    hold on
    ax=gca;
    set(ax,'xScale', 'log')
    set(ax,'YScale', 'log')
    plot(s,Ps)
    xlabel('s')
    ylabel('P(s)')
    xlim([10^(3) 10^(8)])
    ylim([10^(-6) 10^(-1)])

    saveas(gcf,append('Ps_bminutum_chr',num2str(i)),'png')

    i
end


function [s Ps] = contact_probability_HiC(cont_count_mat,resolution)
Ps=zeros(size(cont_count_mat,1)-1,1);
s=1:1:size(cont_count_mat,1)-1;
s=s';
s=s.*resolution;
cont_count=zeros(1,size(cont_count_mat,1)-1);
for i=1:1:size(cont_count_mat,1)-1
    cont_count(i) = sum(diag(cont_count_mat,i));
    Ps(i)=cont_count(i)/(size(cont_count_mat,1)-i);
end
end

