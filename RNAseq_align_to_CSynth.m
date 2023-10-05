clc
clear
get_model = importdata('symbiodinium_microadriaticum_chr2.xyz');
get_RNAseq = importdata('s_microadriaticum_chr2_pilon_strandsum_TPM.bed');

resolution = get_model(2,1)-get_model(1,1);

RNA = [get_RNAseq.data(:,1) get_RNAseq.data(:,2) get_RNAseq.data(:,8)];

[coeff,score] = pca(get_model(:,2:4));

%%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%% use this get_RNAseq.data(:,8)!!! does TPM change as a functino of
%%% radius? or z-axis
%%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

model_PCA = [get_model(:,1) score];

%get xyz coordinates for every base pair by linearly interpolating between the model

model_PCA_interp = zeros(get_model(end,1),4);
model_PCA_interp(:,1) = linspace(1,get_model(end,1),get_model(end,1))';
for i = 1:1:size(model_PCA,1)-1
        model_PCA_interp(resolution*(i-1)+1:resolution*i,2:4) = [linspace(model_PCA(i,2),model_PCA(i+1,2),resolution)' linspace(model_PCA(i,3),model_PCA(i+1,3),resolution)' linspace(model_PCA(i,4),model_PCA(i+1,4),resolution)'];
end

gene_loc = [];
gene_centre_coord = [];

%store coordinates for each base pair in each gene
for i = 1:1:size(RNA,1)
    gene_loc = [gene_loc; model_PCA_interp(RNA(i,1):RNA(i,2),2) model_PCA_interp(RNA(i,1):RNA(i,2),3) model_PCA_interp(RNA(i,1):RNA(i,2),4)];
    gene_centre_coord = [gene_centre_coord; mean(model_PCA_interp(RNA(i,1):RNA(i,2),2)) mean(model_PCA_interp(RNA(i,1):RNA(i,2),3)) mean(model_PCA_interp(RNA(i,1):RNA(i,2),4))];
    size(RNA,1)-i
end

[theta_TPM,rho_TPM] = cart2pol(gene_centre_coord(:,2),gene_centre_coord(:,3));

figure
scatter(rho_TPM, abs(get_RNAseq.data(:,8)))

figure
scatter(gene_centre_coord(:,1), abs(get_RNAseq.data(:,8)))

figure
hold on
histogram(model_PCA_interp(:,2), 25,'Normalization','pdf')
histogram(gene_loc(:,1), 25,'Normalization','pdf')
hold off

[theta,rho] = cart2pol(gene_loc(:,2),gene_loc(:,3));
[theta_model_PCA_interp,rho_model_PCA_interp] = cart2pol(model_PCA_interp(:,3),model_PCA_interp(:,4));

figure
hold on
histogram(rho_model_PCA_interp,25,'Normalization','pdf')
histogram(rho,25,'Normalization','pdf')
legend()
hold off

figure
scatter(theta,rho)

figure
hold on
plot3(model_PCA_interp(:,2),model_PCA_interp(:,3),model_PCA_interp(:,4),'Color', 'g')
plot3(gene_loc(:,1),gene_loc(:,2),gene_loc(:,3),'Color', 'r')
plot3(model_PCA(:,2),model_PCA(:,3),model_PCA(:,4),'Color', 'k')
hold off

numPoints = size(get_model,1);
MyColor = linspace(1,numPoints,numPoints)';
% create a connectivity matrix
Faces = [1:(numPoints-1); 2:numPoints]';

f=figure
hold on
%plot3(get_model(:,2),get_model(:,3),get_model(:,4),'Color', [.6 .6 .6])
%plot3(score(:,1),score(:,2),score(:,3),'Color', [.6 .6 .6])


colormap jet
axis equal
xlim([min(get_model(:,2))*2 max(get_model(:,2))*2])
ylim([min(get_model(:,3))*2 max(get_model(:,3))*2])
zlim([min(get_model(:,4))*2 max(get_model(:,4))*2])
caxis([min(MyColor) max(MyColor)])
c = colorbar;
c.Label.String = 'primary sequence';
patch('Faces', Faces,'Vertices', get_model(:,2:4) ,'FaceColor', 'none', 'FaceVertexCData', MyColor ,'EdgeColor','interp' ,'LineWidth',5, 'FaceAlpha',.5,'EdgeAlpha',.5);
view(30,30)

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
if nargin < 1, m = size(get(gcf,'colormap'),1); end
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b];
end

