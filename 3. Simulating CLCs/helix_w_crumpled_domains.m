% Simulate Helix with Crumpled Domains
clc; clear;

% Parameters
n_domains = 10; % Number of crumpled domains
monomers_per_domain = 400; % domain size (primary sequence)
radius = 4; %radius
turns = 3; % Number of helical turns
helix_height = 20;
domain_size_range = [1.5, 2]; % domain size parameter (3D)
resolution = 5000; %5 kb per monomer

t = linspace(0, 2*pi*turns, n_domains);

% Helix coordinates
helix_x = radius .* cos(t);
helix_y = radius .* sin(t);
helix_z = linspace(0,helix_height,n_domains);

% Accumulate all points
all_points = [];

for i = 1:n_domains
    % Center of current domain
    center = [helix_x(i), helix_y(i), helix_z(i)];

    % Domain size
    r = unifrnd(domain_size_range(1), domain_size_range(2));

    % Number of points in domain
    n_points = round(unifrnd(monomers_per_domain*0.7, monomers_per_domain*1.3));

    % Crumpled blob using Gaussian noise + distortion
    domain = randn(n_points, 3);
    domain = domain .* (0.5 + rand(1,3));  % Irregular distortion
    domain = domain * r;

    % position domain center on helix backbone
    domain = domain + center;

    % save coordinates
    all_points = [all_points; domain];
end


% Plot
% Colormap based on order
n_points = size(all_points, 1);
colors = hsv(n_points); % or 'jet', 'turbo', etc.

% Plot
figure;
scatter3(all_points(:,1), all_points(:,2), all_points(:,3), 200, colors, "filled");
hold on

t = linspace(0, 2*pi*turns, n_points);

% Helix coordinates
helix_x = radius .* cos(t);
helix_y = radius .* sin(t);
helix_z = linspace(0,helix_height,n_points);
% Plot thick interpolated helix path
hold on;
plot3(helix_x*1.5, helix_y*1.5, helix_z, 'k', 'LineWidth', 20);
axis equal;
xlabel('X', 'FontSize', 24);
ylabel('Y', 'FontSize', 24);
zlabel('Z', 'FontSize', 24);
set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
grid off;

D = pdist(all_points); %in microns

D = squareform(D);

%HiC contact probablities from distance
P = (1)./D.^4;
P(isinf(P)) = 1; %self contact probabilities are 1
P(P>1) = 1; %no contact probabilities above 1

figure
imagesc(P);
hold on;
colorbar
xlabel('Primary sequence [bp]', 'fontsize', 24)
ylabel('Primary sequence [bp]', 'fontsize', 24)
xlim([0 size(P,1)])
ylim([0 size(P,1)])
ax = gca;
ax.XTickLabel = ax.XTick*resolution;
ax.YTickLabel = ax.YTick*resolution;
set(gca,'ColorScale','log')
maj_axis_chr=colorbar;
maj_axis_chr.Label.String = 'Contact Probability';
maj_axis_chr.FontSize = 18;
caxis([10^(-4) 10^0]);