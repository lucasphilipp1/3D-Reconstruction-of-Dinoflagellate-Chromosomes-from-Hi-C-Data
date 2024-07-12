clc
clear

num_chroms = 99;

Pm=importdata('s_microadriaticum_chr1_pilon.txt');
avg_diag_Pm = [];
avg_diag_Pm = [avg_diag_Pm; nanmean(diag(Pm,0))];
avg_diag_Pm = [avg_diag_Pm; nanmean(diag(Pm,1))];

Pm_temp = Pm;
Pm_temp_2 = Pm;

Pm_temp(idiag(size(Pm_temp),-1)) = NaN(size(Pm_temp,1)-1,1);
Pm_temp(idiag(size(Pm_temp),0)) = NaN(size(Pm_temp,1),1);
Pm_temp(idiag(size(Pm_temp),1)) = NaN(size(Pm_temp,1)-1,1);
avg_diag_Pm = [avg_diag_Pm; nanmean(nanmean(Pm_temp))];

for i = 1:1:num_chroms
    PCSynth =[];
    Pk=importdata(sprintf('s_kawagutii_HiC_scaffold_%i.txt',i));
    avg_diag_Pk = [];
    avg_diag_Pk = [avg_diag_Pk; nanmean(diag(Pk,0))];
    avg_diag_Pk = [avg_diag_Pk; nanmean(diag(Pk,1))];
    Pk_temp = Pk;
    Pk_temp_2 = Pk;

    Pk_temp(idiag(size(Pk_temp),-1)) = NaN(size(Pk_temp,1)-1,1);
    Pk_temp(idiag(size(Pk_temp),0)) = NaN(size(Pk_temp,1),1);
    Pk_temp(idiag(size(Pk_temp),1)) = NaN(size(Pk_temp,1)-1,1);
    avg_diag_Pk = [avg_diag_Pk; nanmean(nanmean(Pk_temp))];

    Pk = Pk./avg_diag_Pk(3).*avg_diag_Pm(3);

    Pk(idiag(size(Pk),0)) = diag(Pk_temp_2,0)./avg_diag_Pk(1).*avg_diag_Pm(1);
    Pk(idiag(size(Pk),-1)) = diag(Pk_temp_2,-1)./avg_diag_Pk(2).*avg_diag_Pm(2);
    Pk(idiag(size(Pk),1)) = diag(Pk_temp_2,1)./avg_diag_Pk(2).*avg_diag_Pm(2);

    %output for CSynth
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PCSynth = zeros(nat_sum(size(Pk,1)),3);
    count = 0;
    for k=1:1:size(Pk,1)
        for j=k:1:size(Pk,1)
            count = count + 1;
            %PCSynth(count,1)=k;
            %PCSynth(count,2)=j;
            PCSynth(count,1)=k*5000;
            PCSynth(count,2)=j*5000;
            PCSynth(count,3)=Pk(k,j);
        end
    end

    PCSynth(:,3) = round(PCSynth(:,3),3,"significant");
    writematrix(PCSynth,sprintf('s_kawagutii_V3_HiC_scaffold_%i_for_CSynth.txt',i),'Delimiter','tab')
    i
end

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