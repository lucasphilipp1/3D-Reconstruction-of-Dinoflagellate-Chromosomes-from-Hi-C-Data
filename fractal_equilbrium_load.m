clc
clear

%f = fullfile('myfolder','mysubfolder','myfile.m')
source_dir = ('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/GSE18199_Globules/equilibrium/');
N_equil=size(dir([source_dir, '/*.dat']),1);
file_suffix = cell(N_equil,1);

equilibrium_struct=tdfread('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/GSE18199_Globules/equilibrium/equilibrium1.dat');
equilibrium_names = fieldnames(equilibrium_struct);
steps_equil = size(getfield(equilibrium_struct,equilibrium_names{1}),1);
equilibrium = zeros(steps_equil,3,N_equil);

for i=1:1:N_equil
    file_suffix{i} = append(append('equilibrium',num2str(i)),'.dat');
    equilibrium_struct=tdfread(append('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/GSE18199_Globules/equilibrium/',file_suffix{i}));
    equilibrium_names = fieldnames(equilibrium_struct);
    equilibrium(:,:,i) = getfield(equilibrium_struct,equilibrium_names{1});
    i
end

source_dir = ('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/GSE18199_Globules/fractal/');
N_frac=size(dir([source_dir, '/*.dat']),1);
file_suffix = cell(N_frac,1);

fractal_struct=tdfread('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/GSE18199_Globules/fractal/fractal1.dat');
fractal_names = fieldnames(fractal_struct);
steps_frac = size(getfield(fractal_struct,fractal_names{1}),1);
fractal = zeros(steps_frac,3,N_frac);

for i=1:1:N_frac
    file_suffix{i} = append(append('fractal',num2str(i)),'.dat');
    fractal_struct=tdfread(append('/Users/lucasphilipp/Desktop/Research/GitHub/Dinoflagellate_Cholesteric_Large_Files_zenodo/GSE18199_Globules/fractal/',file_suffix{i}));
    fractal_names = fieldnames(fractal_struct);
    fractal(:,:,i) = getfield(fractal_struct,fractal_names{1});
    i
end

N_equil=100;
N_frac=100;

Ps_equil=zeros(steps_equil-1,N_equil);
Ps_equil_average = zeros(steps_equil-1,1);

s_equil=1:1:steps_equil-1;
s_equil=s_equil';

%decide on a distance threshold
num=0;
for k = 1:steps_equil-1
num=num+norm(equilibrium(k,:,1)-equilibrium(k+1,:,1));
end
num=num/(steps_equil-1);
distance_threshold = num*3;

for k = 1:N_equil
    [s_equil Ps_equil(:,k)]=contact_probability_xyz(equilibrium(:,:,k),distance_threshold);
    k
end

for i=1:1:steps_equil-1
    num=0;
    for k = 1:N_equil
        num = num + Ps_equil(i,k);
    end
    Ps_equil_average(i) = num/N_equil;
end

Ps_frac=zeros(steps_frac-1,N_frac);
Ps_frac_average = zeros(steps_frac-1,1);

s_frac=1:1:steps_frac-1;
s_frac=s_frac';

%decide on a distance threshold
num=0;
for k = 1:steps_frac-1
num=num+norm(fractal(k,:,1)-fractal(k+1,:,1));
end
num=num/(steps_frac-1);
distance_threshold = num*3;

for k = 1:N_frac
    [s_frac Ps_frac(:,k)]=contact_probability_xyz(fractal(:,:,k),distance_threshold);
    k
end

for i=1:1:steps_frac-1
    num=0;
    for k = 1:N_frac
        num = num + Ps_frac(i,k);
    end
    Ps_frac_average(i) = num/N_frac;
end

figure
hold on
ax=gca;
set(ax,'xScale', 'log')
set(ax,'YScale', 'log')

A=7.5;
B=4;

xA=linspace(3,100,100);
xB=linspace(5,1000,100);
plot(xA,A*xA.^(-3/2)) %equilibrium globule theory
plot(xB,B*xB.^(-1)) %fractal globule theory

plot(s_equil,Ps_equil_average)
plot(s_frac,Ps_frac_average)

xlabel('s')
ylabel('P(s)')

legend({"-3/2","-1","equilibrium globule","fractal globule"}, 'Location', 'southwest')

%plot fractal globules and equilibrium globules

% figure
% for i = 1:1:1000
% plot3(fractal(:,1,i),fractal(:,2,i),fractal(:,3,i))
% pause(0.1)
% end

