%Author: Lucas Philipp
clc
clear

%data from:
%Zhang, B., & Wolynes, P. G. (2016). Shape transitions and chiral symmetry breaking in the energy landscape of the mitotic chromosome. Physical review letters, 116(24), 248101.

matrix = load('contact_probability_normalized_by_individual.mat');

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

%code to plot eigenvalues and eigenvectors of HiC matrix
[V,D] = eig(matrix.cp);
plot(1:1:size(V(:,1)),V(:,1))
plot(1:1:size(V(:,1)),V(:,2))

%eigenvalues are in D

function sum = nat_sum(x)
sum=0;
for i=1:x
    sum=sum+i;
end
end

