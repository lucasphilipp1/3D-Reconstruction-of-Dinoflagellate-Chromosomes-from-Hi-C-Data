clc
clear

%author: Lucas Philipp
%CSynth requires tab delimited .txt files
%for some reason this bash script:
%for i in $(seq 1 94);
%do cooler dump -t pixels -r chr${i} -r2 chr${i} -b --no-balance -o symbiodinium_microadriaticum_coccoid_chr${i}_temp.txt --join GSE152150_Dino-HiC-cD1plus-R1.smic1.0.mapq_30.1000.mcool::/resolutions/1000; cut -f 2,5,7 symbiodinium_microadriaticum_coccoid_chr${i}_temp.txt | column -ts $'\t' > symbiodinium_microadriaticum_coccoid_chr${i}.txt; rm -f symbiodinium_microadriaticum_coccoid_chr${i}_temp.txt; done
%
%for i in $(seq 1 94);
%do cooler dump -t pixels -r chr${i} -r2 chr${i} -b --no-balance -o symbiodinium_microadriaticum_mastigote_chr${i}_temp.txt --join GSE152150_Dino-HiC-mD1plus-R1.smic1.0.mapq_30.1000.mcool::/resolutions/1000; cut -f 2,5,7 symbiodinium_microadriaticum_mastigote_chr${i}_temp.txt | column -ts $'\t' > symbiodinium_microadriaticum_mastigote_chr${i}.txt; rm -f symbiodinium_microadriaticum_mastigote_chr${i}_temp.txt; done
%outputs space separated .txt files that CSynth cannot read
%this matlab script rewrites these files, making them tab delimited
%USE THIS CODE IMMEADIATELY AFTER ABOVE TERMINAL COMMANDS

%coccoid
contacts = tdfread('/Users/lucasphilipp/Downloads/cholesteric_CSynth_large.txt');
savepath = '/Users/lucasphilipp/Downloads/cholesteric_MATLAB.txt';
matrix=cell2mat(struct2cell(contacts));
writematrix(matrix,savepath,'Delimiter','tab')
size(unique(matrix(:,1)),1)
%16384 is the maximum number of monomers supported by CSynth

%coccoid
% for i=1:94
%     chr_num=num2str(i);
%     filepath = append('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/hic contact probabilities microadriaticum/coccoid/symbiodinium_microadriaticum_coccoid_chr',chr_num);
%     filepath = append(filepath,'.txt');
%     contacts = tdfread(filepath);
%     savepath = append('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/hic contact probabilities microadriaticum/coccoid/symbiodinium_microadriaticum_coccoid_chr',chr_num);
%     savepath = append(savepath,'.txt');
%     matrix=cell2mat(struct2cell(contacts));
%     writematrix(matrix,savepath,'Delimiter','tab')
%     size(unique(matrix(:,1)),1)
%     %16384 is the maximum number of monomers supported by CSynth
% end

% %mastigote
% for i=1:94
%     chr_num=num2str(i);
%     filepath = append('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/hic contact probabilities microadriaticum/mastigote/symbiodinium_microadriaticum_mastigote_chr',chr_num);
%     filepath = append(filepath,'.txt');
%     contacts = tdfread(filepath);
%     savepath = append('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/hic contact probabilities microadriaticum/mastigote/symbiodinium_microadriaticum_mastigote_chr',chr_num);
%     savepath = append(savepath,'.txt');
%     matrix=cell2mat(struct2cell(contacts));
%     writematrix(matrix,savepath,'Delimiter','tab')
%     size(unique(matrix(:,1)),1)
%     %16384 is the maximum number of monomers supported by CSynth
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc
% clear
% contacts = tdfread('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/Positive Control/IMR90/5kb_resolution_intrachromosomal/chr20/MAPQG0/chr20_5kb.RAWobserved');
% names = fieldnames(contacts);
% matrix=zeros(size(getfield(contacts,names{1}),1),3);
% matrix(:,1) = getfield(contacts,names{1});
% matrix(:,2) = getfield(contacts,names{2});
% matrix(:,3) = getfield(contacts,names{3});
% writematrix(matrix,'/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/Positive Control/IMR90/5kb_resolution_intrachromosomal/chr20/MAPQG0/chr20_5kb.RAWobserved.txt','Delimiter','tab')
% size(unique(matrix(:,1)),1)
% %16384 is the maximum number of monomers supported by CSynth
%
% contacts = tdfread('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/Positive Control/IMR90/5kb_resolution_intrachromosomal/chr21/MAPQG0/chr21_5kb.RAWobserved');
% names = fieldnames(contacts);
% matrix=zeros(size(getfield(contacts,names{1}),1),3);
% matrix(:,1) = getfield(contacts,names{1});
% matrix(:,2) = getfield(contacts,names{2});
% matrix(:,3) = getfield(contacts,names{3});
% writematrix(matrix,'/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/Positive Control/IMR90/5kb_resolution_intrachromosomal/chr21/MAPQG0/chr21_5kb.RAWobserved.txt','Delimiter','tab')
% size(unique(matrix(:,1)),1)
% %16384 is the maximum number of monomers supported by CSynth
%
% contacts = tdfread('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/Positive Control/IMR90/5kb_resolution_intrachromosomal/chr22/MAPQG0/chr22_5kb.RAWobserved');
% names = fieldnames(contacts);
% matrix=zeros(size(getfield(contacts,names{1}),1),3);
% matrix(:,1) = getfield(contacts,names{1});
% matrix(:,2) = getfield(contacts,names{2});
% matrix(:,3) = getfield(contacts,names{3});
% writematrix(matrix,'/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric/Positive Control/IMR90/5kb_resolution_intrachromosomal/chr22/MAPQG0/chr22_5kb.RAWobserved.txt','Delimiter','tab')
% size(unique(matrix(:,1)),1)
% %16384 is the maximum number of monomers supported by CSynth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



