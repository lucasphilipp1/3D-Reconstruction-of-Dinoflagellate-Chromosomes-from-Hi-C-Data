%Author: Lucas Philipp
clc
clear

%data from:
%Zhang, B., & Wolynes, P. G. (2016). Shape transitions and chiral symmetry breaking in the energy landscape of the mitotic chromosome. Physical review letters, 116(24), 248101.

matrix = load('contact_probability_smoothed.mat');

%output for CSynth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCSynth = zeros(nat_sum(size(matrix.cp,1)),3);
count = 0;
for i=1:1:size(matrix.cp,1)
    for j=i:1:size(matrix.cp,1)
        count = count + 1;
        PCSynth(count,1)=i*10000;
        PCSynth(count,2)=j*10000;
        PCSynth(count,3)=round(matrix.cp(i,j),3,"significant");
    end
end

PCSynth(PCSynth < 0) = 0;
writematrix(PCSynth,'Naumova_2013_human_metaphase_CSynth.txt','Delimiter','tab')

%smoothing seems to dramatically affect the eigenvalues
%code to plot eigenvalues and eigenvectors of HiC matrix
%eigenvalues are in D

[V,D] = eig(matrix.cp);
figure
plot(1:1:size(V(:,1)),V(:,1))

figure
imagesc(matrix.cp);
colorbar
xlabel('~5kbp monomer index', 'fontsize', 24) %fix this!!!
ylabel('~5kbp monomer index', 'fontsize', 24)
%fix ticks!!!
set(gca,'ColorScale','log')
maj_axis_chr=colorbar;
maj_axis_chr.Label.String = 'Monomer Contact Probability';
maj_axis_chr.FontSize = 18;

matrix = load('contact_probability_normalized_by_individual.mat');

figure
imagesc(matrix.cp);
colorbar
xlabel('~5kbp monomer index', 'fontsize', 24) %fix this!!!
ylabel('~5kbp monomer index', 'fontsize', 24)
%fix ticks!!!
set(gca,'ColorScale','log')
maj_axis_chr=colorbar;
maj_axis_chr.Label.String = 'Monomer Contact Probability';
maj_axis_chr.FontSize = 18;

[V,D] = eig(matrix.cp);
figure
plot(1:1:size(V(:,1)),V(:,1))

function sum = nat_sum(x)
sum=0;
for i=1:x
    sum=sum+i;
end
end

