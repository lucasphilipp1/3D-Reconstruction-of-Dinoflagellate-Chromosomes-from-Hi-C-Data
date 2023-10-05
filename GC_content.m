%requires MATLAB's Bioinformatics Toolbox
clc
clear

% [header_Smic11N,sequence_Smic11N] = fastaread('GSE152150_Smic1.1N.fa');
% Density_Smic11N = cell(size(header_Smic11N,2),1);
% NTStruct_Smic11N = cell(size(header_Smic11N,2),1);
% unmappable_Smic11N = cell(size(header_Smic11N,2),1);
% percent_N_Smic11N = zeros(size(header_Smic11N,2),1);
% 
% % for i = 1:size(header_Smic11N,2)
% %     Density_Smic11N{i} = ntdensity(sequence_Smic11N{i},'Window',5000);
% % end
% 
% for i = 1:size(header_Smic11N,2)
%     NTStruct_Smic11N{i} = basecount(sequence_Smic11N{i},'Ambiguous', 'individual');
%     percent_N_Smic11N(i) = NTStruct_Smic11N{i}.N./length(sequence_Smic11N{i})*100;
% end
% 
% count = logspace(log10(50),log10(100*1000),15);
% for j = 1:1:size(header_Smic11N,2)
%     %seqwordcount does not count overlapping patterns multiple times
%     %however, Ns are considered an any of: A, T, C, G
%     %searching for 'NNN' using seqwordcount counts AAA, NNN, etc...
%     %A, T, C, G are changed to an unrecognized letter so consecutive Ns are counted explicitly
%     %recognized letters: A, C, G, T, R, Y, K, M, S, W, B, D, H, V, N
%     temp_seq=sequence_Smic11N{j};
%     temp_seq = strrep(temp_seq,'A','Z');
%     temp_seq = strrep(temp_seq,'T','Z');
%     temp_seq = strrep(temp_seq,'C','Z');
%     temp_seq = strrep(temp_seq,'G','Z');
% 
%     nhist=[];
%     for i = 1:1:size(count,2)
%         word = repmat('N',1,round(count(i)));
%         nhist=[nhist; seqwordcount(temp_seq, word)];
%     end
% 
%     for i = 1:1:size(count,2)
%         nhist(i)=nhist(i).*round(count(i))./length(sequence_Smic11N{j}).*100;
%     end
% 
%     unmappable_Smic11N{j} = nhist;
%     j
% end
% 
% %percentage of chromosome unmappable given average read-length
% figure
% hold on
% for i = 1:1:size(header_Smic11N,2)
%     plot(count,unmappable_Smic11N{i}', 'Color', [0 0.4470 0.7410])
% end
% %ChIP-seq read length: ~100 bp
% %DiMeLo-seq read length: ~100 000 bp
% xline(100,'-r',{'Short read sequencing','e.g. ChIP-seq & ATAC-seq'}, 'FontSize', 12);
% xline(20*1000,'-g',{'Long read sequencing','e.g. DiMeLo-seq & Fiber-seq'}, 'FontSize', 12);
% ylabel('% chromosome unmappable', 'FontSize', 18)
% xlabel('average readlength [bp]', 'FontSize', 18)
% xlim([0 20*10^4])
% ylim([0 10])
% set(gca, 'XScale', 'log')
% %title('Smic1.1N Symbiodinium microadriaticum', 'FontSize', 18)
% title('Symbiodinium microadriaticum', 'FontSize', 18)

% [header_Smic10,sequence_Smic10] = fastaread('GSE152150_Smic1.0.fa');
% Density_Smic10 = cell(size(header_Smic10,2),1);
% for i = 1:size(header_Smic10,2)
%     Density_Smic10{i} = ntdensity(sequence_Smic10{i},'Window',5000);
% end

[header_SSB01,sequence_SSB01] = fastaread('SSB01_HiC_assembly.fa');
Density_SSB01 = cell(size(header_SSB01,2),1);
NTStruct_SSB01 = cell(size(header_SSB01,2),1);
unmappable_SSB01 = cell(size(header_SSB01,2),1);
percent_N_SSB01 = zeros(size(header_SSB01,2),1);

% for i = 1:size(header_SSB01,2)
%     Density_SSB01{i} = ntdensity(sequence_SSB01{i},'Window',5000);
% end

for i = 1:size(header_SSB01,2)
    NTStruct_SSB01{i} = basecount(sequence_SSB01{i},'Ambiguous', 'individual');
    percent_N_SSB01(i) = NTStruct_SSB01{i}.N./length(sequence_SSB01{i})*100;
end

count = logspace(log10(50),log10(100*1000),15);
for j = 1:1:100
    %seqwordcount does not count overlapping patterns multiple times
    %however, Ns are considered an any of: A, T, C, G
    %searching for 'NNN' using seqwordcount counts AAA, NNN, etc...
    %A, T, C, G are changed to an unrecognized letter so consecutive Ns are counted explicitly
    %recognized letters: A, C, G, T, R, Y, K, M, S, W, B, D, H, V, N
    temp_seq=sequence_SSB01{j};
    temp_seq = strrep(temp_seq,'A','Z');
    temp_seq = strrep(temp_seq,'T','Z');
    temp_seq = strrep(temp_seq,'C','Z');
    temp_seq = strrep(temp_seq,'G','Z');

    nhist=[];
    for i = 1:1:size(count,2)
        word = repmat('N',1,round(count(i)));
        nhist=[nhist; seqwordcount(temp_seq, word)];
    end

    for i = 1:1:size(count,2)
        nhist(i)=nhist(i).*round(count(i))./length(sequence_SSB01{j}).*100;
    end

    unmappable_SSB01{j} = nhist;
    j
end

%percentage of genome unmappable using average read-length
figure
hold on
for i = 1:100
    plot(count,unmappable_SSB01{i}', 'Color', [0 0.4470 0.7410])
end
%ChIP-seq read length: ~100 bp
%DiMeLo-seq read length: ~100 000 bp
xline(100,'-r',{'Short read sequencing','e.g. ChIP-seq & ATAC-seq'}, 'FontSize', 12);
xline(20*1000,'-g',{'Long read sequencing','e.g. DiMeLo-seq & Fiber-seq'}, 'FontSize', 12);
ylabel('% chromosome unmappable', 'FontSize', 18)
xlabel('average readlength [bp]', 'FontSize', 18)
xlim([0 20*10^4])
ylim([0 10])
set(gca, 'XScale', 'log')
%title('Breviolum minutum', 'FontSize', 18)



