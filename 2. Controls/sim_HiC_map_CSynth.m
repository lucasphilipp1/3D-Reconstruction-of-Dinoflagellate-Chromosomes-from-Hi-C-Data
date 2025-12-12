%Fig S13 A & G in paper

chromosome = importdata('Su_Cell_2020_chr21_CSynth_output_SP_0_CF_60_PP_-4.xyz');
input_HiC = importdata('Hi-C_contacts_chromosome21.csv');

resolution = 50000; %bp per bin/monomer

%Hi-C & FISH data is not for the entire chromosome, only about ~75% of the total chromosome
%CSynth expects first monomer to start at 1bp
rows = input_HiC(:,1);
rows(1) = [];
input_HiC(:,1) = [];
input_HiC(1,:) = [];
sto=10400001:50000:46650001;
remove=setxor(sto,rows);
remove=(remove-10400001)./resolution;

figure
imagesc(input_HiC./max(max(input_HiC)));
hold on;
colorbar
xlabel('primary sequence [bp]', 'fontsize', 24)
ylabel('primary sequence [bp]', 'fontsize', 24)
ax = gca;
ax.XTickLabel = ax.XTick*resolution;
ax.YTickLabel = ax.YTick*resolution;
set(gca,'ColorScale','log')
maj_axis_chr=colorbar;
maj_axis_chr.Label.String = 'Contact Probability';
maj_axis_chr.FontSize = 18;
caxis([10^-4 10^0]);

%simulate Hi-C contact map for human IMR90 chr 21
P = 1./squareform(pdist(chromosome(:,2:4))).^4;
P=P./nanmedian(nanmedian(P));

P(:,remove+1) = [];
P(remove+1,:) = [];

figure
imagesc(P./10^4);
hold on;
colorbar
xlabel('primary sequence [bp]', 'fontsize', 24)
ylabel('primary sequence [bp]', 'fontsize', 24)
ax = gca;
ax.XTickLabel = ax.XTick*resolution;
ax.YTickLabel = ax.YTick*resolution;
set(gca,'ColorScale','log')
maj_axis_chr=colorbar;
maj_axis_chr.Label.String = 'Contact Probability';
maj_axis_chr.FontSize = 18;
caxis([10^-4 10^0]);

%omit main and 1st off diagonal for relative error calculation
input_HiC(idiag(size(input_HiC),-1)) = NaN(size(input_HiC,1)-1,1);
input_HiC(idiag(size(input_HiC),0)) = NaN(size(input_HiC,1),1);
input_HiC(idiag(size(input_HiC),1)) = NaN(size(input_HiC,1)-1,1);

P(idiag(size(P),-1)) = NaN(size(P,1)-1,1);
P(idiag(size(P),0)) = NaN(size(P,1),1);
P(idiag(size(P),1)) = NaN(size(P,1)-1,1);

input_HiC(isnan(input_HiC))=0;
P(isnan(P))=0;

input_HiC(isinf(input_HiC))=0;
P(isinf(P))=0;

%relative error
mean(mean((input_HiC-P)./P,'omitnan'),'omitnan')*100

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