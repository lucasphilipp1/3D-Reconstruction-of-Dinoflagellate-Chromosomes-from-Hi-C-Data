clc
clear
%Author: Lucas Philipp
%REQUIRES STATISTICS AND MACHINE LEARNING TOOLBOX!

%%%% USE CELL ARRAY TO STORE MATRICIES/VECTORS OF DIFFERENT SIZES

%%%% TO DO: layer stitch = clostest monomer connection between layers if no
%%%% extra-chromosomal loops
%%%% b (longitudinal component) monomer order
%%%% # OF MONOMERS TAKES ACCOUNTS FOR EXTRACHROMOSOMAL LOOPLENGTH
%%%% loops have average contour length \lambda, gaps between loops have an
%%%% average contour length g, both exponentially distributed

model=[]; %position of DNA fibres
chol_layers = 20;
num_mon = 150; %target for number of monomers in the thicket layer

%ellipse parameters for overall chromosome profile
min_axis_chr=0.5; %in microns

spacing = min_axis_chr*sqrt(pi/num_mon); %mesh spacing in microns
%2D grid of points
x = -1:spacing:1; %mesh in microns
y = -1:spacing:1; %mesh in microns

maj_axis_chr=chol_layers*spacing; %in microns

%ellipse minor/major axes ratio for a single layer, for circle set ratio=1
%major axis of layer ellipse = minor axis of chromosome ellipse
ell_ratio = 0.90;

%longitudinal component magnitude
b = 0.130 * 0.025; %in microns

%ensure ~num_mon monomers in the thickest layer (where r=b)
%7.655Mbp/(10 layers)∗1monomer/5kbp≈(150 monomers)/layer
mon_per_layer=zeros(chol_layers,1);

%cholesteric pitch P: length in microns along long axis of the chromosome
%corresponding to a full turn of the nematic director
pitch = 1.30*2; %in microns

dtheta_layer=360*spacing/pitch; %DNA fibres in a layer are rotated dtheta_layer degrees counterclockwise relative to the previous layer

%for vector field describing orientation of DNA fibres
lin_vec_x =[];
lin_vec_y =[];

intra_disc_loop_potential_starts = [];

i=1;
for z = 0-(spacing*(chol_layers-1))/2:spacing:spacing*(chol_layers-1)-(spacing*(chol_layers-1))/2 %thickest/middle part of chromosome is z=0

    [X,Y] = meshgrid(x,y);
    r = (min_axis_chr/maj_axis_chr)*sqrt(maj_axis_chr^2-z^2); %major axis of disk at this z position

    %layer mask
    indicator = sqrt((X./r).^2 + (Y./(r.*ell_ratio)).^2) - 1 < 0;
    xEllipse = X(indicator);
    yEllipse = Y(indicator);

    %every second line of monomers y->-y so ordering does a "weave"
    [c,ia,ic] = unique(xEllipse);
    flip_these=[];
    for k = 1:2:length(c)
        flip_these = [flip_these; find(ic==k)];
    end
    yEllipse(flip_these) = -1.*yEllipse(flip_these);

    R=[cosd((i-1)*dtheta_layer) -sind((i-1)*dtheta_layer); sind((i-1)*dtheta_layer) cosd((i-1)*dtheta_layer)]; %create 2D rotation matix, applies to (X,Y) pairs
    rotcoord=[];
    for j=1:1:size(xEllipse,1)
        rotcoord=[rotcoord; [xEllipse(j), yEllipse(j)]*R']; %multiply by rotation matrix
    end

    intra_disc_loop_potential_starts = [intra_disc_loop_potential_starts; size(model,1) + flip_these];

    model = [model; [rotcoord(:,1), rotcoord(:,2), z.*ones(size(xEllipse))]]; %append positions of DNA fibres for new layer

    lin_vec_x = [lin_vec_x; spacing*cosd((i-1)*dtheta_layer).*ones(size(xEllipse))]; %append orientation of DNA fibres for new layer
    lin_vec_y = [lin_vec_y; spacing*sind((i-1)*dtheta_layer).*ones(size(xEllipse))]; %append orientation of DNA fibres for new layer

    mon_per_layer(i)=size(xEllipse,1);
    i=i+1;
end

%orientation of DNA fibres, no longituindal component: fvec=(cos(az), sin(az), 0)
lin_vec_0=ones(size(model));
%orientation of DNA fibres, longitudnal component=b fvec=(cos(az), sin(az), b) %normalize to unit length???
lin_vec_b=ones(size(model));

lin_vec_0(:,1) = lin_vec_x;
lin_vec_0(:,2) = lin_vec_y;
lin_vec_0(:,3) = zeros(size(model(:,3)));

lin_vec_b(:,1) = lin_vec_x;
lin_vec_b(:,2) = lin_vec_y;
lin_vec_b(:,3) = b.*ones(size(model(:,3))); %not affected by normalization

%another way to taper: no taper, elliptical, hockey rink
z = 0-(spacing*(chol_layers-1))/2:spacing:spacing*(chol_layers-1)-(spacing*(chol_layers-1))/2;
r = (min_axis_chr/maj_axis_chr)*sqrt(maj_axis_chr^2-z.^2);

%DNA FIBRES HAVE NO LONGITUDINAL COMPONENT
% figure
% set(gcf, 'Position', get(0, 'Screensize'));
% hold on
% q0=quiver3(model(:,1)-lin_vec_0(:,1)./2,model(:,2)-lin_vec_0(:,2)./2,model(:,3)-lin_vec_0(:,3)./2,lin_vec_0(:,1),lin_vec_0(:,2),lin_vec_0(:,3),'off'); %orientation vector rotates about midpoint instead of base
% q0.ShowArrowHead = 'off';
% %plot3(model(:,1),model(:,2),model(:,3),'.r');
% plot3(r./sqrt(2),r./sqrt(2),z,'.k')
% %[L,M,N] = cylinder(r);
% %N=N.*chol_layer_sep*chol_layers;
% %N=N-(chol_layer_sep*(chol_layers-1))/2;
% %surf(L,M,N)
% %xlim([-1.5 1.5])
% %ylim([-1.5 1.5])
% xlim([-chol_layers*chol_layer_sep chol_layers*chol_layer_sep])
% ylim([-chol_layers*chol_layer_sep chol_layers*chol_layer_sep])
% zlim([-chol_layers*chol_layer_sep chol_layers*chol_layer_sep])
% xlabel('x','FontSize', 24)
% ylabel('y','FontSize', 24)
% zlabel('z','FontSize', 24)
% title('Cholesteric Model of a Dinoflagellate Chromosome', 'FontSize', 24)
% axis vis3d
% view(0,0)
% for i=0:1:360-1
%     view(i,0)
%     pause(0.01);
% end
%
% view(0,10)

%DNA FIBRES HAVE LONGITUDINAL COMPONENT
% figure
% set(gcf, 'Position', get(0, 'Screensize'));
% hold on
% qb=quiver3(model(:,1)-lin_vec_b(:,1)./2,model(:,2)-lin_vec_b(:,2)./2,model(:,3)-lin_vec_b(:,3)./2,lin_vec_b(:,1),lin_vec_b(:,2),lin_vec_b(:,3),'off'); %orientation vector rotates about midpoint instead of base
% qb.ShowArrowHead = 'off';
% %plot3(model(:,1),model(:,2),model(:,3),'.r');
% plot3(r./sqrt(2),r./sqrt(2),z,'.k')
% %xlim([-1.5 1.5])
% %ylim([-1.5 1.5])
% xlim([-chol_layers*chol_layer_sep chol_layers*chol_layer_sep])
% ylim([-chol_layers*chol_layer_sep chol_layers*chol_layer_sep])
% zlim([-chol_layers*chol_layer_sep chol_layers*chol_layer_sep])
% xlabel('x','FontSize', 24)
% ylabel('y','FontSize', 24)
% zlabel('z','FontSize', 24)
% %title('Cholesteric Model of a Dinoflagellate Chromosome', 'FontSize', 24)
% axis vis3d
% for i=0:1:2*360-1
%     view(i,10)
%     pause(0.01);
% end
%
% %rotate entire model and vector field by phi
% phi=45; %degrees
% Rphi = [1 0 0; 0 cosd(phi) -sind(phi); 0 sind(phi) cosd(phi)];
% model_rot=model*Rphi;
% lin_vec_0_rot=lin_vec_0*Rphi;
% lin_vec_b_rot=lin_vec_b*Rphi;
% figure
% set(gcf, 'Position', get(0, 'Screensize'));
% q0=quiver3(model_rot(:,1)-lin_vec_0_rot(:,1)./2,model_rot(:,2)-lin_vec_0_rot(:,2)./2,model_rot(:,3)-lin_vec_0_rot(:,3)./2,lin_vec_0_rot(:,1),lin_vec_0_rot(:,2),lin_vec_0_rot(:,3),'off'); %orientation vector rotates about midpoint instead of base
% q0.ShowArrowHead = 'off';
%
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% zlim([-0.15 0.15])
% view(0,90)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% title('Rotated model for oblique cross-section viewing')

%INTER DISK LOOPS
idx=find(diff(model(:,3))>0); %new disc starts when z coordinate changes
inter_disc_loop_length = 10; %requested loop length along each axis
%looplength = looplength*3;

%make sure the accepted loop length is close to the requested loop length
num_inter_disc_loops = size(idx,1);
inter_disc_loops = zeros(inter_disc_loop_length, 3, num_inter_disc_loops);
disc_starts = zeros(num_inter_disc_loops+1, 3);
disc_ends = zeros(num_inter_disc_loops+1, 3);
save_disc_starts = zeros(num_inter_disc_loops+1, 3); %remains off-latice, matches model coordinates
save_disc_ends = zeros(num_inter_disc_loops+1, 3); %remains off-latice, matches model coordinates

for i=1:1:num_inter_disc_loops+1
    if i == 1
        disc_ends(i,:)=model(idx(1),:);
        disc_starts(i,:)=model(1,:);
    elseif i == num_inter_disc_loops+1
        disc_ends(i,:)=model(end,:);
        disc_starts(i,:)=model(idx(end)+1,:);
    else
        disc_ends(i,:)=model(idx(i),:);
        disc_starts(i,:)=model(idx(i-1)+1,:);
    end
end

for i=1:1:num_inter_disc_loops+1
    save_disc_starts(i,:)=disc_starts(i,:);
    save_disc_ends(i,:)=disc_ends(i,:);

    disc_starts(i,:)=spacing*round(disc_starts(i,:)./spacing);
    disc_ends(i,:)=spacing*round(disc_ends(i,:)./spacing);
end

count = 1;
for i=1:1:num_inter_disc_loops
    %adjacent disks connected
    if mod(count,2)==1
        inter_disc_loops(:,:,i)=constrained_self_avoiding_RW_3D(disc_starts(count,:),disc_starts(count+1,:),inter_disc_loop_length,spacing);
    else
        inter_disc_loops(:,:,i)=constrained_self_avoiding_RW_3D(disc_ends(count,:),disc_ends(count+1,:),inter_disc_loop_length,spacing);
    end
    for j=1:1:inter_disc_loop_length
        if sqrt((inter_disc_loops(j,1,i)).^2 + (inter_disc_loops(j,2,i)).^2) <= min_axis_chr
            [theta,rho] = cart2pol(inter_disc_loops(j,1,i),inter_disc_loops(j,2,i));
            rho=rho+abs(2*(min_axis_chr-sqrt((inter_disc_loops(j,1,i)).^2 + (inter_disc_loops(j,2,i)).^2)));
            [inter_disc_loops(j,1,i),inter_disc_loops(j,2,i)] = pol2cart(theta,rho);
        end
    end
    %resolved problem where two intertior parts of the loop are being reflected to opposite ends of the cylinder, the loop then crosses through the cylinder
    %reject walks with any steps of length spacing*3 or greater
    j=1;
    while j<=size(inter_disc_loops(:,:,i),1)-1
        while sqrt((inter_disc_loops(j,1,i)-inter_disc_loops(j+1,1,i)).^2 + (inter_disc_loops(j,2,i)-inter_disc_loops(j+1,2,i)).^2 + (inter_disc_loops(j,3,i)-inter_disc_loops(j+1,3,i)).^2) > spacing*2
            if mod(count,2)==1
                inter_disc_loops(:,:,i)=constrained_self_avoiding_RW_3D(disc_starts(count,:),disc_starts(count+1,:),inter_disc_loop_length,spacing);
                j=1;
            else
                inter_disc_loops(:,:,i)=constrained_self_avoiding_RW_3D(disc_ends(count,:),disc_ends(count+1,:),inter_disc_loop_length,spacing);
                j=1;
            end
            for k=1:1:inter_disc_loop_length
                if sqrt((inter_disc_loops(k,1,i)).^2 + (inter_disc_loops(k,2,i)).^2) <= min_axis_chr
                    [theta,rho] = cart2pol(inter_disc_loops(k,1,i),inter_disc_loops(k,2,i));
                    rho=rho+abs(2*(min_axis_chr-sqrt((inter_disc_loops(k,1,i)).^2 + (inter_disc_loops(k,2,i)).^2)));
                    [inter_disc_loops(k,1,i),inter_disc_loops(k,2,i)] = pol2cart(theta,rho);
                end
            end
        end
        j=j+1;
    end
    count = count + 1;
end

[Xcyl,Ycyl,Zcyl] = cylinder(min_axis_chr);
Zcyl = Zcyl*(max(model(:,3))-min(model(:,3)));
Zcyl = Zcyl - max(model(:,3));

%%% CHANGE NAME: NUMLOOPS, LOOPS 

% f=figure
% count=1;
% for i=1:1:numloops
%     clf(f)
%     hold on
%     surf(Xcyl,Ycyl,Zcyl,'FaceAlpha',0.2);
% 
%     if mod(count,2)==1
%         plot3(disc_starts(count,1),disc_starts(count,2),disc_starts(count,3),'o','Color','g')
%         plot3(disc_starts(count+1,1),disc_starts(count+1,2),disc_starts(count+1,3),'o','Color','r')
%         count = count + 1;
%     else
%         plot3(disc_ends(count,1),disc_ends(count,2),disc_ends(count,3),'o','Color','g')
%         plot3(disc_ends(count+1,1),disc_ends(count+1,2),disc_ends(count+1,3),'o','Color','r')
%         count = count + 1;
%     end
% 
%     plot3(loops(:,1,i),loops(:,2,i),loops(:,3,i),'-','Color','k')
%     plot3(loops(1,1,i),loops(1,2,i),loops(1,3,i),'o','MarkerFaceColor','g')
%     plot3(loops(end,1,i),loops(end,2,i),loops(end,3,i),'o','MarkerFaceColor','r')
% 
%     xlim([min(model(:,1))*2 max(model(:,1))*2])
%     ylim([min(model(:,2))*2 max(model(:,2))*2])
%     zlim([min(model(:,3))*2 max(model(:,3))*2])
%     view(0,90)
%     pause(0.1)
% end
% hold off

% figure
% hold on
% count=1;
% for i=1:1:num_inter_disc_loops
%     if mod(count,2)==1
%         plot3(disc_starts(count,1),disc_starts(count,2),disc_starts(count,3),'o','Color','g')
%         plot3(disc_starts(count+1,1),disc_starts(count+1,2),disc_starts(count+1,3),'o','Color','r')
%         count = count + 1;
%     else
%         plot3(disc_ends(count,1),disc_ends(count,2),disc_ends(count,3),'o','Color','g')
%         plot3(disc_ends(count+1,1),disc_ends(count+1,2),disc_ends(count+1,3),'o','Color','r')
%         count = count + 1;
%     end
% 
%     plot3(loops(:,1,i),loops(:,2,i),loops(:,3,i),'-','Color','k')
%     plot3(loops(1,1,i),loops(1,2,i),loops(1,3,i),'o','MarkerFaceColor','g')
%     plot3(loops(end,1,i),loops(end,2,i),loops(end,3,i),'o','MarkerFaceColor','r')
% 
%     xlim([min(model(:,1))*2 max(model(:,1))*2])
%     ylim([min(model(:,2))*2 max(model(:,2))*2])
%     zlim([min(model(:,3))*2 max(model(:,3))*2])
%     view(30,30)
% end
% plot3(model(:,1),model(:,2),model(:,3))

% insert loops into primary sequence in model
idx=zeros(num_inter_disc_loops,1);
idx=find(diff(model(:,3))>0); %new disc starts when z coordinate changes

model_unflipped=model;

model(1:idx(1),:) = flipud(model(1:idx(1),:));
for i=1:1:size(idx,1)
    if mod(i,2)==0
        if i == size(idx,1)
            model(idx(i)+1:end,:) = flipud(model(idx(i)+1:end,:));
        else
            model(idx(i)+1:idx(i+1),:) = flipud(model(idx(i)+1:idx(i+1),:));
        end
    end
end

intra_disc_loop_potential_starts = sort(intra_disc_loop_potential_starts);

append=[];
for i = 1:1:size(intra_disc_loop_potential_starts,1)
    [match, idx2] = ismember(model_unflipped(intra_disc_loop_potential_starts(i),:), model, 'rows');
    intra_disc_loop_potential_starts(i)=idx2;
    append = [append; idx2+1];
    append = [append; idx2-1];
end
intra_disc_loop_potential_starts = [intra_disc_loop_potential_starts; append];

intra_disc_loop_potential_starts = unique(intra_disc_loop_potential_starts,'rows');

intra_disc_loop_potential_starts = sort(intra_disc_loop_potential_starts);

intra_disc_loop_potential_starts(intra_disc_loop_potential_starts == 0, :) = [];
intra_disc_loop_potential_starts(intra_disc_loop_potential_starts == size(model,1)+1, :) = [];

save=[];
for i = 3:1:size(intra_disc_loop_potential_starts,1)-2
    if ismember(model(intra_disc_loop_potential_starts(i)-2,:), model(intra_disc_loop_potential_starts,:), 'rows') & ismember(model(intra_disc_loop_potential_starts(i)+2,:), model(intra_disc_loop_potential_starts,:), 'rows')
        save=[save; i];
    end
end

intra_disc_loop_potential_starts(save)=[];

% Preallocate output
new_model = zeros(size(model,1)+inter_disc_loop_length*num_inter_disc_loops,3);
% Find indices for old data
addRows = ismember(1:size(model,1), idx)*inter_disc_loop_length;

oldDataInd = (1:size(model,1)) + cumsum([0, addRows(1:end-1)]);

% Add in old data
new_model(oldDataInd,:) = model(:,:);

%%% i like this way of adding in new data better
% temp=find(all(new_new_model==0,2));
% for i=1:1:num_intra_disc_loops
%     % Find indices for new data
%     newDataInd = temp((i-1)*intra_disc_loop_length+1:i*intra_disc_loop_length);
%     % Add in new data
%     new_new_model(newDataInd,:) = intra_disc_loops(:,:,i);
% end

for i=1:1:num_inter_disc_loops
    % Find indices for new data
    newDataInd = idx(i)+1+inter_disc_loop_length*(i-1):1:idx(i)+inter_disc_loop_length*i;
    % Add in new data
    new_model(newDataInd,:) = inter_disc_loops(:,:,i);
end

for i = 1:1:size(intra_disc_loop_potential_starts,1)
    [match, idx3] = ismember(model(intra_disc_loop_potential_starts(i),:), new_model, 'rows');
    intra_disc_loop_potential_starts(i)=idx3;
end

%INTRA DISK LOOPS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

intra_disc_loop_length = 4;
num_intra_disc_loops = 20;

%intra_disc_loop_potential_starts are indecies for new_model
%idx4 can be removed. potentially idx5

%choose 20 random discs
new_model(intra_disc_loop_potential_starts(i),3)

%choose positions within those discs

%make sure no overlap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx=find(diff(model(:,3))>0); %new disc starts when z coordinate changes
intra_disc_loops = zeros(intra_disc_loop_length, 3, num_intra_disc_loops);

halfrange = size(intra_disc_loop_potential_starts,1)/ 2;
idx4 = (randperm(halfrange, num_intra_disc_loops) - 1) * 2 + 1;

intra_disc_loop_starts=new_model(intra_disc_loop_potential_starts(idx4),:);
intra_disc_loop_ends=new_model(intra_disc_loop_potential_starts(idx4+1),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:1:num_intra_disc_loops
    intra_disc_loop_starts(i,:)=spacing*round(intra_disc_loop_starts(i,:)./spacing);
    intra_disc_loop_ends(i,:)=spacing*round(intra_disc_loop_ends(i,:)./spacing);
end

for i=1:1:num_intra_disc_loops
    intra_disc_loops(:,:,i)=constrained_self_avoiding_RW_3D(intra_disc_loop_starts(i,:),intra_disc_loop_ends(i,:),intra_disc_loop_length,spacing);
    for j=1:1:intra_disc_loop_length
        if sqrt((intra_disc_loops(j,1,i)).^2 + (intra_disc_loops(j,2,i)).^2) <= min_axis_chr
            [theta,rho] = cart2pol(intra_disc_loops(j,1,i),intra_disc_loops(j,2,i));
            rho=rho+abs(2*(min_axis_chr-sqrt((intra_disc_loops(j,1,i)).^2 + (intra_disc_loops(j,2,i)).^2)));
            [intra_disc_loops(j,1,i),intra_disc_loops(j,2,i)] = pol2cart(theta,rho);
        end
    end
    %resolved problem where two intertior parts of the loop are being reflected to opposite ends of the cylinder, the loop then crosses through the cylinder
    %reject walks with any steps of length spacing*3 or greater
    j=1;
    while j<=size(intra_disc_loops(:,:,i),1)-1
        while sqrt((intra_disc_loops(j,1,i)-intra_disc_loops(j+1,1,i)).^2 + (intra_disc_loops(j,2,i)-intra_disc_loops(j+1,2,i)).^2 + (intra_disc_loops(j,3,i)-intra_disc_loops(j+1,3,i)).^2) > spacing*2
            intra_disc_loops(:,:,i)=constrained_self_avoiding_RW_3D(intra_disc_loop_starts(i,:),intra_disc_loop_ends(i,:),intra_disc_loop_length,spacing);
            j=1;
            for k=1:1:intra_disc_loop_length
                if sqrt((intra_disc_loops(k,1,i)).^2 + (intra_disc_loops(k,2,i)).^2) <= min_axis_chr
                    [theta,rho] = cart2pol(intra_disc_loops(k,1,i),intra_disc_loops(k,2,i));
                    rho=rho+abs(2*(min_axis_chr-sqrt((intra_disc_loops(k,1,i)).^2 + (intra_disc_loops(k,2,i)).^2)));
                    [intra_disc_loops(k,1,i),intra_disc_loops(k,2,i)] = pol2cart(theta,rho);
                end
            end
        end
        j=j+1;
    end
    i
end

idx5=intra_disc_loop_potential_starts(idx4);

% Preallocate output
new_new_model = zeros(size(new_model,1)+intra_disc_loop_length*num_intra_disc_loops,3);
% Find indices for old data
addRows = ismember(1:size(new_model,1), intra_disc_loop_potential_starts(idx4))*intra_disc_loop_length;

oldDataInd = (1:size(new_model,1)) + cumsum([0, addRows(1:end-1)]);

% Add in old data
new_new_model(oldDataInd,:) = new_model(:,:);

temp=find(all(new_new_model==0,2));
for i=1:1:num_intra_disc_loops
    % Find indices for new data
    newDataInd = temp((i-1)*intra_disc_loop_length+1:i*intra_disc_loop_length);
    % Add in new data
    new_new_model(newDataInd,:) = intra_disc_loops(:,:,i);
end

figure
hold on
plot3(new_new_model(:,1),new_new_model(:,2),new_new_model(:,3))
for i = 1:1:num_intra_disc_loops
plot3(intra_disc_loops(:,1,i),intra_disc_loops(:,2,i),intra_disc_loops(:,3,i),'o','MarkerFaceColor','g')
end
xlim([min(new_new_model(:,1))*2 max(new_new_model(:,1))*2])
ylim([min(new_new_model(:,2))*2 max(new_new_model(:,2))*2])
zlim([min(new_new_model(:,3))*2 max(new_new_model(:,3))*2])
xlabel('x','FontSize', 24)
ylabel('y','FontSize', 24)
zlabel('z','FontSize', 24)
view(30,30)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numPoints = size(new_new_model,1);

MyColor = linspace(1,numPoints,numPoints)';
% create a connectivity matrix
Faces = [1:(numPoints-1); 2:numPoints]';

% skip=5;
% 
% f=figure
% hold on
% plot3(new_new_model(:,1),new_new_model(:,2),new_new_model(:,3),'Color', [.6 .6 .6])
% colormap jet
% axis equal
% xlim([min(new_new_model(:,1))*2 max(new_new_model(:,1))*2])
% ylim([min(new_new_model(:,2))*2 max(new_new_model(:,2))*2])
% zlim([min(new_new_model(:,3))*2 max(new_new_model(:,3))*2])
% caxis([min(MyColor) max(MyColor)])
% c = colorbar;
% c.Label.String = 'primary sequence';
% for i=skip+1:skip:numPoints
%     clf(f)
%     plot3(new_new_model(:,1),new_new_model(:,2),new_new_model(:,3),'Color', [.6 .6 .6])
%     xlim([min(new_new_model(:,1))*2 max(new_new_model(:,1))*2])
%     ylim([min(new_new_model(:,2))*2 max(new_new_model(:,2))*2])
%     zlim([min(new_new_model(:,3))*2 max(new_new_model(:,3))*2])
%     view(30,30)
%     caxis([min(MyColor) max(MyColor)])
%     c = colorbar;
%     c.Label.String = 'primary sequence';
%     patch('Faces', Faces(i-skip:i-1,:) ,'Vertices', new_new_model(1:i,:) ,'FaceColor', 'none', 'FaceVertexCData', MyColor(1:i,:) ,'EdgeColor','interp' ,'LineWidth',5, 'FaceAlpha',.5,'EdgeAlpha',.5);
%     pause(0.1);
% end

%%%
%y=y_0 + 1/(a sin(phi)*cos(phi))*log(abs(sin(x*a*sin(phi))+btan(phi))
%%%

%euclidian distance between monomers
%D = pdist(model); %in microns
D = pdist(new_new_model); %in microns

D = squareform(D);

%HiC contact probablities from distance
P = spacing./D;
P(isinf(P)) = 1; %self contact probabilities are 1

%are these coming from disc/loop connections?
P(P>1) = 1; %no contact probabilities above 1

%output for GEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%writematrix(P,'cholesteric_GEM_HiC.txt','Delimiter','tab')

P_loci = zeros(size(P,1),1);
count = 0;
for i=1:1:size(P,1)
    count = count + 1;
    P_loci(count,1)=i*5000;
end

%writematrix(P_loci,'cholesteric_GEM_loci.txt','Delimiter','tab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output for CSynth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCSynth = zeros(nat_sum(size(P,1)),3);
count = 0;
for i=1:1:size(P,1)
    for j=i:1:size(P,1)
        count = count + 1;
        PCSynth(count,1)=i;
        PCSynth(count,2)=j;
        PCSynth(count,3)=P(i,j);
    end
end

%writematrix(PCSynth,'cholesteric_CSynth.txt','Delimiter','tab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %plot hic matrix
% figure
% imagesc(P);
% colorbar
% xlabel('~5kbp monomer index', 'fontsize', 24) %fix this!!!
% ylabel('~5kbp monomer index', 'fontsize', 24)
% %fix ticks!!!
% maj_axis_chr=colorbar;
% maj_axis_chr.Label.String = 'Monomer Contact Probability';
% maj_axis_chr.FontSize = 18;
% %what are the boxes in this plot???

%extra-chromosomal loops
%add variable amounts of sequence between: lines in a disc, between discs,
%jumble disc ordering

%local interactions should be interpretable even if configuration of extrachromosomal
%looping is unknown

%what about bigger and smaller circles? not ellipses?

%%%GENERATE AN ENSEBLE OF CHOLESTERIC CONFORMATIONS WITH DIFFERENT
%%%EXTRACHROMSOMAL LOOPS
%%%REAL HIC DATA IS AN ENSEMBLE OF CHROMOSOME CONFORMATIONS OF MANY CELLS

%total number of monomers
tot_monomers_in_body = sum(mon_per_layer);

%(The possible number of pairs of genomic positions separated by d on a given chromosome is Lc-d, where Lc is the length of the chromosome.)
distance_threshold = spacing*5;
%[s Ps] = contact_probability_xyz(model,distance_threshold);
[s Ps] = contact_probability_xyz(new_model,distance_threshold);

figure
hold on
ax=gca;
xlim([0,400]);
plot(s,Ps)
legendStrings = "d_{cutoff} = " + string(distance_threshold) + "um";
legend(legendStrings, 'Location', 'southwest')
xlabel('s in multiples of 5kbp')
ylabel('P(s)')
set(ax,'xScale', 'log')
set(ax,'YScale', 'log')
grid on

function sum = nat_sum(x)
sum=0;
for i=1:x
    sum=sum+i;
end
end

%
% %P(s) calculation.
% %P(s) plots were calculated in two ways.
% %First, P(s) was calculated genome-wide using valid chromatin interaction pairs
% %with the following script from the cooltools package with all default settings (https://github.com/mirnylab/ cooltools).
% %For P(s) calculations for single chromosomes at the level of chromosomal Hi-C domains,
% %we used Hi-C data binned and balanced at 50-kb resolution. Hi-C domains borders were calculated by insulation analysis (see above).
% %The grid of domain borders defines a set of squares throughout the Hi-C interaction map. P(s) was calculated for each square by
% %plotting the average of each diagonal of 50-kb bins in the square as a function of s.
% %For all squares not centered at the main diagonal, the values of P(s) for the smallest and largest s were left out because
% %they are calculated only for 1 bin and thus are noisy.
%
% %does a CLC model with a particular parameter set give rise to a unique scaling exponent?
% %or can one scaling exponent correspond to multiple geometries?

