get_model = importdata('human_metaphase.xyz');
get_model(:,1,:) = [];

number_of_monomers = size(get_model,1);



%To further smooth out noise, we determined the position of each genomic locus as the center
%of mass of the nearest in sequence 4 Mb genomic loci, i.e.
%(i − 2Mb) : (i + 2Mb), around that locus.

%keeps track of the summation to give the H value for the entire structure
total_H = 0;

%the resolution of the human metaphase HiC data is 10000 bp

%use a sequence separation of 25*10^6 bp
T = number_of_monomers/8;

%where E and F are the midpoints of the vectors −−→AB
%and −−→CD respectively

M = zeros(8,2);
count=1;

%do the minus version: get_model(i-round(5/4*T),:,:);

%iterate through the entire cholesteric model
for i=1:1:round(5/4*T)
    %get monomers
    point_A = get_model(i,:,:);
    point_B = get_model(i+round(1/2*T),:,:);
    point_C = get_model(i+round(3/4*T),:,:);
    point_D = get_model(i+round(5/4*T),:,:);

    %get vectors
    point_AB = point_B-point_A; %B-A
    point_CD = point_D-point_C; %D-C
    point_EF = point_CD-point_AB; %CD-AB

    %calculate H of ABCD
    calculate_H = (dot(point_EF, cross(point_CD, point_AB)))/(norm(point_AB)*norm(point_CD)*norm(point_EF));
    
    % if isnan(calculate_H)
    %     point_A
    %     point_B
    %     point_C
    %     point_D
    %     point_AB
    %     point_CD
    %     point_EF
    %     calculate_H
    % end

    %if NaN is encountered, ignore
    if isnan(calculate_H)
        calculate_H = 0;
    end

    total_H = total_H + calculate_H;
    if rem(count,100)==0
        M(count/100,:) = [count total_H];
    end
    count = count+1;
end
plot(M(:,1),M(:,2));

%H > 0 therefore right hand helix


