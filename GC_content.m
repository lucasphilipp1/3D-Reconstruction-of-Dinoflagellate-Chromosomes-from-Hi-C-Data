%requires MATLAB's Bioinformatics Toolbox
clc
clear

[header_Smic11N,sequence_Smic11N] = fastaread('GSE152150_Smic1.1N.fa');

Density_Smic11N = cell(size(header_Smic11N,2),1);
for i = 1:size(header_Smic11N,2)
    Density_Smic11N{i} = ntdensity(sequence_Smic11N{i},'Window',5000);
end

[header_Smic10,sequence_Smic10] = fastaread('GSE152150_Smic1.0.fa');
Density_Smic10 = cell(size(header_Smic10,2),1);
for i = 1:size(header_Smic10,2)
    Density_Smic10{i} = ntdensity(sequence_Smic10{i},'Window',5000);
end

[header_SSB01,sequence_SSB01] = fastaread('SSB01_HiC_assembly.fa');
Density_SSB01 = cell(size(header_SSB01,2),1);
for i = 1:size(header_SSB01,2)
    Density_SSB01{i} = ntdensity(sequence_SSB01{i},'Window',5000);
end