% Loads IGM polymer coordinates and plots the 3D chromosome structure.
% Color the polymer by primary sequence position or by Topological Associating Domain coordinates.
%
% Usage:
%   1. Set file path to IGM coordinate file.
%   2. Set COLOR_MODE to 'sequence' or 'tad'.
%   3. Run the script.

clc; clear;

% --- File paths ---
COORDS_FILE = 'path_to_IGM_coordinate_file';   % IGM polymer coordinates (.txt)

% --- Genome parameters ---
RESOLUTION = 5000;               % Hi-C resolution in bp/monomer
CHROMOSOME_LENGTH_BP = 19282064; % Total chromosome length in bp
                                 % Common values:
                                 %   s_microadriaticum chr1: 19282064
                                 %   f_kawagutii chr1:       17515114
                                 %   b_minutum chr1:         11576343

% --- Coloring mode ---
% 'sequence' : color by position along primary sequence (jet colormap)
% 'tad'      : color each TAD domain a distinct color
COLOR_MODE = 'sequence';

% --- TAD boundaries in base pairs [start, stop] ---
% Add or remove rows as needed.
TAD_BP = [
      21449    3136855;   % TAD 1
    3138645    5456633;   % TAD 2
    5440909    7815782;   % TAD 3
    7817867    8660723;   % TAD 4
    8677664   10342262;   % TAD 5
   10355695   14413730;   % TAD 6
   14436643   19268194;   % TAD 7
];

% --- Plot options ---
SPLINE_UPSAMPLING = 3;   % Smoothing factor for spline interpolation (1 = no smoothing)
LINE_WIDTH        = 2;   % Polymer line width

%% ============================================================
chromosome = importdata(COORDS_FILE);
xyz = chromosome(:, 1:3);

% Align to principal axes via PCA
[~, xyz] = pca(xyz);

% Smooth with cubic spline interpolation
N  = size(xyz, 1);
t  = 1:N;
tq = linspace(1, N, N * SPLINE_UPSAMPLING);

xq = spline(t, xyz(:,1), tq);
yq = spline(t, xyz(:,2), tq);
zq = spline(t, xyz(:,3), tq);
xyz_smooth = [xq', yq', zq'];

fig = figure('Color', 'w');
pos = get(fig, 'Position');
set(fig, 'Position', pos .* [1 1 1.4 1]);
hold on; axis equal; axis off;

N_smooth = size(xyz_smooth, 1);

switch lower(COLOR_MODE)

    %% -- Color by primary sequence position --
    case 'sequence'
        cmap = jet(N_smooth - 1);
        for j = 1:N_smooth-1
            plot3(xyz_smooth(j:j+1, 1), ...
                  xyz_smooth(j:j+1, 2), ...
                  xyz_smooth(j:j+1, 3), ...
                  'Color', cmap(j,:), 'LineWidth', LINE_WIDTH);
        end

        colormap(jet);
        clim([0, N - 1]);
        c = colorbar('Location', 'eastoutside');
        c.Label.String = 'Primary Sequence [Mb]';
        c.FontSize = 24;
        c.TickLabels = arrayfun(@(x) sprintf('%.0f', x * RESOLUTION / 1e6), ...
                                c.Ticks, 'UniformOutput', false);

    %% -- Color by TAD --
    case 'tad'
        total_monomers = round(CHROMOSOME_LENGTH_BP / RESOLUTION);
        TAD_MON = round(TAD_BP ./ CHROMOSOME_LENGTH_BP * total_monomers);
        TAD_MON = TAD_MON * SPLINE_UPSAMPLING;  % scale to upsampled indices

        % Color palette (cycles if more TADs than colors)
        color_palette = [
            0.85  0.15  0.15;   % red
            0.15  0.60  0.15;   % green
            0.20  0.40  0.85;   % blue
            0.10  0.70  0.65;   % teal
            0.75  0.15  0.75;   % magenta
            0.90  0.65  0.10;   % orange
            0.50  0.50  0.50;   % grey
        ];

        num_tads   = size(TAD_MON, 1);
        num_colors = size(color_palette, 1);
        leg_handles = gobjects(num_tads, 1);

        for k = 1:num_tads
            idx_start = max(1, TAD_MON(k, 1));
            if k < num_tads
                idx_end = min(N_smooth, TAD_MON(k, 2));
            else
                idx_end = N_smooth;   % last TAD runs to end
            end

            col = color_palette(mod(k-1, num_colors) + 1, :);

            leg_handles(k) = plot3(xyz_smooth(idx_start:idx_end, 1), ...
                                   xyz_smooth(idx_start:idx_end, 2), ...
                                   xyz_smooth(idx_start:idx_end, 3), ...
                                   'Color', col, 'LineWidth', LINE_WIDTH);
        end

        legend_labels = arrayfun(@(x) sprintf('TAD %d', x), 1:num_tads, ...
                                 'UniformOutput', false);
        lgd = legend(leg_handles, legend_labels, 'Location', 'northeast');
        lgd.Box = 'off';
        lgd.FontSize = 16;

    otherwise
        error('Unknown COLOR_MODE "%s". Use ''sequence'' or ''tad''.', COLOR_MODE);
end

margin = 0.05;
for dim = 1:3
    vals = xyz_smooth(:, dim);
    span = max(vals) - min(vals);
    xlims(dim,:) = [min(vals) - margin*span, max(vals) + margin*span];
end
xlim(xlims(1,:)); ylim(xlims(2,:)); zlim(xlims(3,:));

view(0, 90);
set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
