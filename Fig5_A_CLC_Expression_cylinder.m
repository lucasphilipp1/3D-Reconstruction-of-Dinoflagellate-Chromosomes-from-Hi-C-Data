clc
clear
% Define the cylinder
t = 0:pi/10:2*pi;
r = 1 + 0.01*cos(t); %subtle wiggle with z to create many patches along z

% Create cylinders
[cx,cy,cz]  = cylinder(r,100);
cz = cz.*2 - 1;

%shift cylinders
cx          = cx(:,ceil(length(cz)/2):end);   % Remove top half 
cy          = cy(:,ceil(length(cz)/2):end);   % of cylinder
cz          = cz(:,ceil(length(cz)/2):end);

% [cx,cy,cz]  = cylinder(r);
% cz = cz.*2 - 1;

%define plane for cross-section
p1 = [-1.35 0.05 -1.15];
p2 = [-1.35 0.05 1.15];
p3 = [1.35 0.05 1.15];  
p4 = [1.35 0.05 -1.15];

x = [p1(1) p2(1) p3(1) p4(1)];
y = [p1(2) p2(2) p3(2) p4(2)];
z = [p1(3) p2(3) p3(3) p4(3)];

% grey core patch
p1 = [-1 0.06 -1];
p2 = [-1 0.06 1];
p3 = [1 0.06 1];  
p4 = [1 0.06 -1];

x2 = [p1(1) p2(1) p3(1) p4(1)];
y2 = [p1(2) p2(2) p3(2) p4(2)];
z2 = [p1(3) p2(3) p3(3) p4(3)];

%define color for each surface patch
C = zeros(size(cz,1),size(cz,2),3);
for i = 1:1:size(cz,1)
    for j = 1:1:size(cz,2)
        if rand(1)>0.5
        C(i,j,:) = [0 0.4470 0.7410]; %blue patches
        else
        C(i,j,:) = [0.8500 0.1250 0.0380]; %red patches
        end
    end
end

%%% create semi-circular top and bottom lid
[sx_top,sy_top,sz_top] = sphere;
sx_top          = sx_top(ceil(length(sz_top)/2):end,:);   % Remove bottom half 
sy_top          = sy_top(ceil(length(sz_top)/2):end,:);   % of sphere
sz_top          = sz_top(ceil(length(sz_top)/2):end,:);
sz_top = sz_top.*0.1;
sz_top = sz_top + 1;

sx_top          = sx_top(:,1:ceil(length(sz_top)/2));   % Remove right half 
sy_top          = sy_top(:,1:ceil(length(sz_top)/2));   % of hemisphere
sz_top          = sz_top(:,1:ceil(length(sz_top)/2));

[sx_bottom,sy_bottom,sz_bottom] = sphere;
sx_bottom       = sx_bottom(1:ceil(length(sz_bottom)/2),:);   % Remove top half 
sy_bottom       = sy_bottom(1:ceil(length(sz_bottom)/2),:);   % of sphere
sz_bottom       = sz_bottom(1:ceil(length(sz_bottom)/2),:);
sz_bottom = sz_bottom.*0.1;
sz_bottom = sz_bottom - 1;

sx_bottom          = sx_bottom(:,1:ceil(length(sz_bottom)/2));   % Remove right half 
sy_bottom          = sy_bottom(:,1:ceil(length(sz_bottom)/2));   % of hemisphere
sz_bottom          = sz_bottom(:,1:ceil(length(sz_bottom)/2));

%define color for each surface patch
C_top = zeros(size(sz_top,1),size(sz_top,2),3);
for i = 1:1:size(sz_top,1)
    for j = 1:1:size(sz_top,2)
        if rand(1)>0.5
        C_top(i,j,:) = [0 0.4470 0.7410]; %blue patches
        else
        C_top(i,j,:) = [0.8500 0.1250 0.0380]; %red patches
        end
    end
end

C_bottom = zeros(size(sz_bottom,1),size(sz_bottom,2),3);
for i = 1:1:size(sz_bottom,1)
    for j = 1:1:size(sz_bottom,2)
        if rand(1)>0.5
        C_bottom(i,j,:) = [0 0.4470 0.7410]; %blue patches
        else
        C_bottom(i,j,:) = [0.8500 0.1250 0.0380]; %red patches
        end
    end
end

figure;
hold on
fill3(x, y, z, [1 1 1], 'facealpha',0);
fill3(x2, y2, z2, [0 0 0], 'facealpha',0,'edgecolor','none');
surf(cx, cy, cz, C, 'facealpha',.6,'edgecolor','none')
surf(sx_top, sy_top, sz_top, C_top, 'facealpha',.6,'edgecolor','none')
surf(sx_bottom, sy_bottom, sz_bottom, C_bottom, 'facealpha',.6,'edgecolor','none')
grid off
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])
axis equal
view(120,-12.5)
xlim([-1.35 1.35])
ylim([-1.35 1.35])
zlim([-1.35 1.35])