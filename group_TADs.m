%group TADs by asphericity

%0<asphericity<0.2
%0.2<asphericity<0.4
%0.4<asphericity<0.6
%0.6<asphericity<0.8
% 
% all_TADs_microadriaticum = importdata('symbiodinium_microadriaticum_allchroms_TADs_asphericity.txt');
% 
% idx=find(0<all_TADs_microadriaticum.data(:,3) & all_TADs_microadriaticum.data(:,3)<0.2);
% T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
% writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_asphericity_0_0.2.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(0.2<all_TADs_microadriaticum.data(:,3) & all_TADs_microadriaticum.data(:,3)<0.4);
% T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
% writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_asphericity_0.2_0.4.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(0.4<all_TADs_microadriaticum.data(:,3) & all_TADs_microadriaticum.data(:,3)<0.6);
% T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
% writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_asphericity_0.4_0.6.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(0.6<all_TADs_microadriaticum.data(:,3) & all_TADs_microadriaticum.data(:,3)<0.8);
% T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
% writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_asphericity_0.6_0.8.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% %%%
% 
% all_TADs_kawagutii = importdata('s_kawagutii_V3_allchroms_TADs_asphericity.txt');
% 
% idx=find(0<all_TADs_kawagutii.data(:,3) & all_TADs_kawagutii.data(:,3)<0.2);
% T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
% writetable(T,'symbiodinium_kawagutii_allchroms_TADs_asphericity_0_0.2.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(0.2<all_TADs_kawagutii.data(:,3) & all_TADs_kawagutii.data(:,3)<0.4);
% T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
% writetable(T,'symbiodinium_kawagutii_allchroms_TADs_asphericity_0.2_0.4.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(0.4<all_TADs_kawagutii.data(:,3) & all_TADs_kawagutii.data(:,3)<0.6);
% T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
% writetable(T,'symbiodinium_kawagutii_allchroms_TADs_asphericity_0.4_0.6.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(0.6<all_TADs_kawagutii.data(:,3) & all_TADs_kawagutii.data(:,3)<0.8);
% T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
% writetable(T,'symbiodinium_kawagutii_allchroms_TADs_asphericity_0.6_0.8.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% %group TADs by size
% %make histogram
% 
% TAD_size = all_TADs_microadriaticum.data(:,2) - all_TADs_microadriaticum.data(:,1);
% 
% figure
% histogram(TAD_size)
% 
% idx=find(0<TAD_size & TAD_size<1*10^6);
% T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
% writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_size_0_1*10^6.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(1*10^6<TAD_size & TAD_size<2*10^6);
% T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
% writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_size_1*10^6_2*10^6.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(2*10^6<TAD_size & TAD_size<3*10^6);
% T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
% writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_size_2*10^6_3*10^6.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(3*10^6<TAD_size & TAD_size<4*10^6);
% T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
% writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_size_3*10^6_4*10^6.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(4*10^6<TAD_size & TAD_size<5*10^6);
% T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
% writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_size_4*10^6_5*10^6.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% TAD_size = all_TADs_kawagutii.data(:,2) - all_TADs_kawagutii.data(:,1);
% 
% figure
% histogram(TAD_size)
% 
% idx=find(0<TAD_size & TAD_size<1*10^6);
% T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
% writetable(T,'symbiodinium_kawagutii_allchroms_TADs_size_0_1*10^6.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(1*10^6<TAD_size & TAD_size<2*10^6);
% T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
% writetable(T,'symbiodinium_kawagutii_allchroms_TADs_size_1*10^6_2*10^6.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(2*10^6<TAD_size & TAD_size<3*10^6);
% T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
% writetable(T,'symbiodinium_kawagutii_allchroms_TADs_size_2*10^6_3*10^6.txt','Delimiter','\t', 'WriteVariableNames',0)
% 
% idx=find(3*10^6<TAD_size & TAD_size<4*10^6);
% T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
% writetable(T,'symbiodinium_kawagutii_allchroms_TADs_size_3*10^6_4*10^6.txt','Delimiter','\t', 'WriteVariableNames',0)

%group TADs into equal sample size groups based on asphericity percentile

all_TADs_microadriaticum = importdata('symbiodinium_microadriaticum_allchroms_TADs_asphericity.txt');

idx=find(0<all_TADs_microadriaticum.data(:,3) & all_TADs_microadriaticum.data(:,3)<0.2001);
T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_asphericity_first_quartile.txt','Delimiter','\t', 'WriteVariableNames',0)

idx=find(0.2001<all_TADs_microadriaticum.data(:,3) & all_TADs_microadriaticum.data(:,3)<0.3058);
T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_asphericity_second_quartile.txt','Delimiter','\t', 'WriteVariableNames',0)

idx=find(0.3058<all_TADs_microadriaticum.data(:,3) & all_TADs_microadriaticum.data(:,3)<0.4286);
T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_asphericity_third_quartile.txt','Delimiter','\t', 'WriteVariableNames',0)

idx=find(0.4286<all_TADs_microadriaticum.data(:,3) & all_TADs_microadriaticum.data(:,3)<0.8);
T = table(all_TADs_microadriaticum.textdata(idx,:),all_TADs_microadriaticum.data(idx,1:2));
writetable(T,'symbiodinium_microadriaticum_allchroms_TADs_asphericity_fourth_quartile.txt','Delimiter','\t', 'WriteVariableNames',0)

%%%

all_TADs_kawagutii = importdata('s_kawagutii_V3_allchroms_TADs_asphericity.txt');

idx=find(0<all_TADs_kawagutii.data(:,3) & all_TADs_kawagutii.data(:,3)<0.3666);
T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
writetable(T,'symbiodinium_kawagutii_allchroms_TADs_asphericity_first_quartile.txt','Delimiter','\t', 'WriteVariableNames',0)

idx=find(0.3666<all_TADs_kawagutii.data(:,3) & all_TADs_kawagutii.data(:,3)<0.4359);
T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
writetable(T,'symbiodinium_kawagutii_allchroms_TADs_asphericity_second_quartile.txt','Delimiter','\t', 'WriteVariableNames',0)

idx=find(0.4359<all_TADs_kawagutii.data(:,3) & all_TADs_kawagutii.data(:,3)<0.5522);
T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
writetable(T,'symbiodinium_kawagutii_allchroms_TADs_asphericity_third_quartile.txt','Delimiter','\t', 'WriteVariableNames',0)

idx=find(0.5522<all_TADs_kawagutii.data(:,3) & all_TADs_kawagutii.data(:,3)<0.8);
T = table(all_TADs_kawagutii.textdata(idx,:),all_TADs_kawagutii.data(idx,1:2));
writetable(T,'symbiodinium_kawagutii_allchroms_TADs_asphericity_fourth_quartile.txt','Delimiter','\t', 'WriteVariableNames',0)

