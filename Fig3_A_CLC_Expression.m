%code to create 3D visual model of surface-localized gene expression on CLC chromosomes
%using a divergent strand-specific expression colormap

clc
clear
% Define the ellipse
x   = 0;    y   = 0;    z   = 0;
tl  = 1.5;   tw  = 1;   td  = 1;

% Create ellipsoid
[ex,ey,ez]  = ellipsoid(x, y, z, tl, tw, td, 40);
ex          = ex(:,ceil(length(ez)/2):end);   % Remove top half 
ey          = ey(:,ceil(length(ez)/2):end);   % of ellipsoid
ez          = ez(:,ceil(length(ez)/2):end);

%[ex,ey,ez]  = ellipsoid(x, y, z, tl, tw, td,40);

[ex2,ey2,ez2]  = ellipsoid(x, y, z, tl, 0, td, 40);
ex2          = ex2(:,floor(length(ez2)/2):end);   % Remove top half 
ey2          = ey2(:,floor(length(ez2)/2):end);   % of ellipsoid
ez2          = ez2(:,floor(length(ez2)/2):end);

%define plane for cross-section
p1 = [-1.75 0.05 -1.15];
p2 = [-1.75 0.05 1.15];
p3 = [1.75 0.05 1.15];  
p4 = [1.75 0.05 -1.15];

x = [p1(1) p2(1) p3(1) p4(1)];
y = [p1(2) p2(2) p3(2) p4(2)];
z = [p1(3) p2(3) p3(3) p4(3)];

%define color for each surface patch
C = zeros(size(ez,1),size(ez,2),3);
for i = 1:1:size(ez,1)
    for j = 1:1:size(ez,2)
        if rand(1)>0.5
        C(i,j,:) = [0 0.4470 0.7410]; %blue patches
        else
        C(i,j,:) = [0.8500 0.1250 0.0380]; %red patches
        end
    end
end

figure;
hold on
fill3(x, y, z, [1 1 1]);
surf(ex, ey, ez, C, 'facealpha',.6,'edgecolor','none')
surf(ex2, ey2, ez2, 'facecolor','black','facealpha',.2,'edgecolor','none')
grid off
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])
axis equal
view(57.5,25)
xlim([-1.75 1.75])
ylim([-1.75 1.75])
zlim([-1.75 1.75])



    

