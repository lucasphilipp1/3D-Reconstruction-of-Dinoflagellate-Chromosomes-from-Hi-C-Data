clc
clear
%Author: Lucas Philipp
%REQUIRES STATISTICS AND MACHINE LEARNING TOOLBOX!

%%%% TO DO:
%%%% inter-disc extra-chromosomal loops, jumble disc ordering
%%%% what about bigger and smaller circles? not ellipses?
%%%% does a CLC model with a particular parameter set give rise to a unique scaling exponent? or can one scaling exponent correspond to multiple geometries?

%%% GENERATE AN ENSEBLE OF CHOLESTERIC CONFORMATIONS WITH DIFFERENT EXTRACHROMSOMAL LOOPS
%%% REAL HIC DATA IS AN ENSEMBLE OF CHROMOSOME CONFORMATIONS OF MANY CELLS
%%% compute P(s) curve directly from HiC map, using cooltools without
%%% distance cutoff

model=[]; %position of DNA fibres
chol_layers = 20;
num_mon = 50; %target for number of monomers in the thicket layer

%ellipse parameters for overall chromosome profile
min_axis_chr=0.5; %in microns

spacing = min_axis_chr*sqrt(pi/num_mon); %mesh spacing in microns
%2D grid of points
x = -1:spacing:1; %mesh in microns
y = -1:spacing:1; %mesh in microns

maj_axis_chr=chol_layers*spacing; %in microns

%ellipse minor/major axes ratio for a single layer, for circle set ratio=1
%major axis of layer ellipse = minor axis of chromosome ellipse
%ell_ratio = 0.90;
ell_ratio = 1;

%longitudinal component magnitude
b = 0.130 * 0.025; %in microns

%ensure ~num_mon monomers in the thickest layer (where r=b)
%7.655Mbp/(10 layers)∗1monomer/5kbp≈(150 monomers)/layer
mon_per_layer=zeros(chol_layers,1);

%cholesteric pitch P: length in microns along long axis of the chromosome
%corresponding to a full turn of the nematic director
pitch = 1.30/1.5; %in microns

dtheta_layer=360*spacing/pitch; %DNA fibres in a layer are rotated dtheta_layer degrees counterclockwise relative to the previous layer

%for vector field describing orientation of DNA fibres
lin_vec_x =[];
lin_vec_y =[];

intra_disc_loop_potential_starts_ind = [];

i=1;
for z = 0-(spacing*(chol_layers-1))/2:spacing:spacing*(chol_layers-1)-(spacing*(chol_layers-1))/2 %thickest/middle part of chromosome is z=0

    [X,Y] = meshgrid(x,y);
    r = (min_axis_chr/maj_axis_chr)*sqrt(maj_axis_chr^2-z^2); %major axis of disk at this z position

    %layer mask
    indicator = sqrt((X./r).^2 + (Y./(r.*ell_ratio)).^2) - 1 < 0;
    xEllipse = X(indicator);
    yEllipse = Y(indicator);

    %every second line of monomers y->-y so ordering does a "weave"
    [c,ia,ib] = unique(xEllipse);
    flip_these=[];
    for k = 1:2:length(c)
        flip_these = [flip_these; find(ib==k)];
    end
    yEllipse(flip_these) = -1.*yEllipse(flip_these);

    R=[cosd((i-1)*dtheta_layer) -sind((i-1)*dtheta_layer); sind((i-1)*dtheta_layer) cosd((i-1)*dtheta_layer)]; %create 2D rotation matix, applies to (X,Y) pairs
    rotcoord=[];
    for j=1:1:size(xEllipse,1)
        rotcoord=[rotcoord; [xEllipse(j), yEllipse(j)]*R']; %multiply by rotation matrix
    end

    %create array of potential start positions for intra-disc extra chromosomal loops
    intra_disc_loop_potential_starts_ind = [intra_disc_loop_potential_starts_ind; size(model,1) + flip_these];

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
%z = 0-(spacing*(chol_layers-1))/2:spacing:spacing*(chol_layers-1)-(spacing*(chol_layers-1))/2;
%r = (min_axis_chr/maj_axis_chr)*sqrt(maj_axis_chr^2-z.^2);

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
%rotate entire model and vector field by phi
phi=45; %degrees
Rphi = [1 0 0; 0 cosd(phi) -sind(phi); 0 sind(phi) cosd(phi)];
model_rot=model*Rphi;
lin_vec_0_rot=lin_vec_0*Rphi;
lin_vec_b_rot=lin_vec_b*Rphi;
% figure
% set(gcf, 'Position', get(0, 'Screensize'));
% q0=quiver3(model_rot(:,1)-lin_vec_0_rot(:,1)./2,model_rot(:,2)-lin_vec_0_rot(:,2)./2,model_rot(:,3)-lin_vec_0_rot(:,3)./2,lin_vec_0_rot(:,1),lin_vec_0_rot(:,2),lin_vec_0_rot(:,3),'off'); %orientation vector rotates about midpoint instead of base
% q0.ShowArrowHead = 'off';
%
% xlim([-1 1])
% ylim([-1 1])
% zlim([-0.15 0.15])
% view(0,90)
% xlabel('x','FontSize', 24)
% ylabel('y','FontSize', 24)
% zlabel('z','FontSize', 24)
% title('Rotated model for oblique cross-section viewing','FontSize', 24)

%INTER DISK LOOPS
idx=find(diff(model(:,3))>0); %new disc starts when z coordinate changes
mean_looplength = 20; %requested loop length along each axis

model_unflipped=model;

%choose potential intra-disc loop start positions, i.e. perimeter points on discs
intra_disc_loop_potential_starts_ind = sort(intra_disc_loop_potential_starts_ind);

append=[];
%RENAME IDX2
for i = 1:1:size(intra_disc_loop_potential_starts_ind,1)
    [match, idx2] = ismember(model_unflipped(intra_disc_loop_potential_starts_ind(i),:), model, 'rows');
    intra_disc_loop_potential_starts_ind(i)=idx2;
    append = [append; idx2+1];
    append = [append; idx2-1];
end
intra_disc_loop_potential_starts_ind = [intra_disc_loop_potential_starts_ind; append];

intra_disc_loop_potential_starts_ind = unique(intra_disc_loop_potential_starts_ind,'rows');

intra_disc_loop_potential_starts_ind = sort(intra_disc_loop_potential_starts_ind);

intra_disc_loop_potential_starts_ind(intra_disc_loop_potential_starts_ind == 0, :) = [];
intra_disc_loop_potential_starts_ind(intra_disc_loop_potential_starts_ind == size(model,1)+1, :) = [];

save=[];
for i = 3:1:size(intra_disc_loop_potential_starts_ind,1)-2
    if ismember(model(intra_disc_loop_potential_starts_ind(i)-2,:), model(intra_disc_loop_potential_starts_ind,:), 'rows') & ismember(model(intra_disc_loop_potential_starts_ind(i)+2,:), model(intra_disc_loop_potential_starts_ind,:), 'rows')
        save=[save; i];
    end
end

intra_disc_loop_potential_starts_ind(save)=[];

num_inter_disc_loops = size(idx,1);
num_intra_disc_loops = 0;

inter_disc_loops = cell(1, num_inter_disc_loops); %cell array
disc_starts = zeros(num_inter_disc_loops+1, 3); %on lattice points
disc_ends = zeros(num_inter_disc_loops+1, 3);  %on lattice points
save_disc_starts = zeros(num_inter_disc_loops+1, 3); %remains off-latice, matches model coordinates
save_disc_ends = zeros(num_inter_disc_loops+1, 3); %remains off-latice, matches model coordinates

%find starts and ends to each disc
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

%KEEP GENERATING LOOPS UNTIL A CERTAIN SEQUENCE LENGTH IS REACHED.
%OTHERWISE DISCS ARE CONNECTED BY A DIRECT VERTICAL STEP BY DEFAULT (CHECK)
%MOVE SAMPLING EXPONENTIAL DIST TO TOP OF THE CODE

%generate loops with looplength as mean
inter_disc_loop_lengths = exprnd(mean_looplength, 1, num_inter_disc_loops);

%if the generated loop length is rounded to 0, remove it from the list
%inter_disc_loop_lengths(find(round(inter_disc_loop_lengths)==0))=[];
for i=1:1:num_inter_disc_loops
    if round(inter_disc_loop_lengths(i)) == 0
        inter_disc_loop_lengths(i) = 1;
    end
end

count = 1;
for i=1:1:size(inter_disc_loop_lengths,2)
    i
    round(inter_disc_loop_lengths(i))
    %adjacent disks connected
    if mod(count,2)==1
        %convert 2D array to cell
        inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_starts(count,:),disc_starts(count+1,:),round(inter_disc_loop_lengths(i)),spacing));
    else
        %convert 2D array to cell
        inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_ends(count,:),disc_ends(count+1,:),round(inter_disc_loop_lengths(i)),spacing));
    end
    for j=1:1:round(inter_disc_loop_lengths(i))
        if sqrt((inter_disc_loops{i}{j,1}).^2 + (inter_disc_loops{i}{j,2}).^2) <= min_axis_chr
            [theta,rho] = cart2pol(inter_disc_loops{i}{j,1},inter_disc_loops{i}{j,2});
            rho=rho+abs(2*(min_axis_chr-sqrt((inter_disc_loops{i}{j,1}).^2 + (inter_disc_loops{i}{j,2}).^2)));
            [inter_disc_loops{i}{j,1},inter_disc_loops{i}{j,2}] = pol2cart(theta,rho);
        end
    end
    %resolved problem where two intertior parts of the loop are being reflected to opposite ends of the cylinder, the loop then crosses through the cylinder
    %reject walks with any steps of length spacing*2 or greater
    j=1;
    while j<=size(inter_disc_loops{i},1)-1
        while sqrt((inter_disc_loops{i}{j,1})-inter_disc_loops{i}{j+1,1}).^2 + (inter_disc_loops{i}{j,2}-inter_disc_loops{i}{j+1,2}).^2 + (inter_disc_loops{i}{j,3}-inter_disc_loops{i}{j+1,3}).^2 > spacing*2
            if mod(count,2)==1
                inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_starts(count,:),disc_starts(count+1,:),round(inter_disc_loop_lengths(i)),spacing));
                j=1;
            else
                inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_ends(count,:),disc_ends(count+1,:),round(inter_disc_loop_lengths(i)),spacing));
                j=1;
            end
            for k=1:1:round(inter_disc_loop_lengths(i))
                if sqrt((inter_disc_loops{i}{k,1}).^2 + (inter_disc_loops{i}{k,2}).^2) <= min_axis_chr
                    [theta,rho] = cart2pol(inter_disc_loops{i}{k,1},inter_disc_loops{i}{k,2});
                    rho=rho+abs(2*(min_axis_chr-sqrt((inter_disc_loops{i}{k,1}).^2 + (inter_disc_loops{i}{k,2}).^2)));
                    [inter_disc_loops{i}{k,1},inter_disc_loops{i}{k,2}] = pol2cart(theta,rho);
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

% f=figure
% fig=gcf;
% fig.Position(3:4)=[1500,1500];
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
% for i=1:1:numloops
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

%by flipping every second disc the start of the next disc will be directly above the end of the disc below
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

total_length_inter_disc = 0;
for i=1:1:num_inter_disc_loops
    total_length_inter_disc = total_length_inter_disc + round(inter_disc_loop_lengths(i));
end

% insert inter-disc loops into primary sequence in model
%%%RENAME NEW MODEL AND NEW_NEW_MODEL

new_model = insert_loops(num_inter_disc_loops, total_length_inter_disc, model, inter_disc_loops, idx);

% figure
% plot3(new_model(:,1),new_model(:,2),new_model(:,3))
% xlim([min(model(:,1))*2 max(model(:,1))*2])
% ylim([min(model(:,2))*2 max(model(:,2))*2])
% zlim([min(model(:,3))*2 max(model(:,3))*2])
% xlabel('x','FontSize', 24)
% ylabel('y','FontSize', 24)
% zlabel('z','FontSize', 24)
% view(30,30)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RENAME IDX3
%this maps the list of indicies for intra-disc loop starts from model -> new_model
for i = 1:1:size(intra_disc_loop_potential_starts_ind,1)
    [match, idx3] = ismember(model(intra_disc_loop_potential_starts_ind(i),:), new_model, 'rows');
    intra_disc_loop_potential_starts_ind(i)=idx3;
end

%ensure intra-disc loops don't start/end at the same position as an inter-disc loop
for i = 1:1:size(save_disc_starts,1)
    common = ismember(save_disc_starts(i,:),new_model(intra_disc_loop_potential_starts_ind,:),'rows');
    if common == 1
        RowIdx = find(ismember(new_model(intra_disc_loop_potential_starts_ind,:), save_disc_starts(i,:),'rows'));
        intra_disc_loop_potential_starts_ind(RowIdx)=[];
    end
    common = ismember(save_disc_ends(i,:),new_model(intra_disc_loop_potential_starts_ind,:),'rows');
    if common == 1
        RowIdx = find(ismember(new_model(intra_disc_loop_potential_starts_ind,:), save_disc_ends(i,:),'rows'));
        intra_disc_loop_potential_starts_ind(RowIdx)=[];
    end
end

% figure
% hold on
% plot3(new_model(:,1),new_model(:,2),new_model(:,3),'-','Color','k')
% plot3(new_model(intra_disc_loop_potential_starts,1),new_model(intra_disc_loop_potential_starts,2),new_model(intra_disc_loop_potential_starts,3),'o','MarkerFaceColor','r')
% xlim([min(new_model(:,1))*2 max(new_model(:,1))*2])
% ylim([min(new_model(:,2))*2 max(new_model(:,2))*2])
% for i = 0:0.05:1
% zlim([max(new_model(:,3))*(i-0.04) max(new_model(:,3))*(i+0.04)])
% pause(5)
% end
% view(0,90)
% hold off

%choose two neighbouring perimeter points on the same disc
num_intra_disc_loops = 10;

%intra_disc_loop_potential_starts are indicies for new_model
cluster=clusterdata(intra_disc_loop_potential_starts_ind,1);
for i = 1:1:max(cluster)
    if size(find(cluster==i),1) > 2
        intra_disc_loop_potential_starts_ind(find(cluster==i)) = [];
    end
end

intra_disc_loop_potential_starts_ind=intra_disc_loop_potential_starts_ind(find(diff(intra_disc_loop_potential_starts_ind)==1));
p = randperm(size(intra_disc_loop_potential_starts_ind,1),num_intra_disc_loops);
intra_disc_loop_starts_ind = intra_disc_loop_potential_starts_ind(p);
intra_disc_loop_ends_ind = intra_disc_loop_starts_ind + 1;

intra_disc_loop_starts = zeros(num_intra_disc_loops+1, 3); %on lattice points
intra_disc_loop_ends = zeros(num_intra_disc_loops+1, 3);  %on lattice points

for i=1:1:num_intra_disc_loops
    intra_disc_loop_starts(i,:)=spacing*round(new_model(intra_disc_loop_starts_ind(i),:)./spacing);
    intra_disc_loop_ends(i,:)=spacing*round(new_model(intra_disc_loop_ends_ind(i),:)./spacing);
end

% figure
% hold on
% plot3(new_model(:,1),new_model(:,2),new_model(:,3))
% plot3(new_model(intra_disc_loop_starts_ind,1),new_model(intra_disc_loop_starts_ind,2),new_model(intra_disc_loop_starts_ind,3),'o','MarkerFaceColor','r')
% plot3(new_model(intra_disc_loop_ends_ind,1),new_model(intra_disc_loop_ends_ind,2),new_model(intra_disc_loop_ends_ind,3),'o','MarkerFaceColor','g')
% xlim([min(new_model(:,1))*2 max(new_model(:,1))*2])
% ylim([min(new_model(:,2))*2 max(new_model(:,2))*2])
% zlim([min(new_model(:,3))*2 max(new_model(:,3))*2])
% xlabel('x','FontSize', 24)
% ylabel('y','FontSize', 24)
% zlabel('z','FontSize', 24)
% view(30,30)

%KEEP GENERATING LOOPS UNTIL A CERTAIN SEQUENCE LENGTH IS REACHED.
%OTHERWISE DISCS ARE CONNECTED BY A DIRECT VERTICAL STEP BY DEFAULT (CHECK)
%MOVE SAMPLING EXPONENTIAL DIST TO TOP OF THE CODE

%generate loops with looplength as mean
intra_disc_loop_lengths = exprnd(mean_looplength, 1, num_intra_disc_loops);

%if the generated loop length is rounded to 0, remove it from the list
%intra_disc_loop_lengths(find(round(intra_disc_loop_lengths)==0))=[];
%!!!!!!!!!!!!!!!!
for i=1:1:num_intra_disc_loops
    if round(intra_disc_loop_lengths(i)) == 0
        intra_disc_loop_lengths(i) = 1;
    end
end

%Generate intra-disc loops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intra_disc_loops = cell(1, num_intra_disc_loops); %cell array

for i=1:1:size(intra_disc_loop_lengths,2)
    intra_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(intra_disc_loop_starts(i,:),intra_disc_loop_ends(i,:),round(intra_disc_loop_lengths(i)),spacing));
    for j=1:1:round(intra_disc_loop_lengths(i))
        if sqrt((intra_disc_loops{i}{j,1}).^2 + (intra_disc_loops{i}{j,2}).^2) <= min_axis_chr
            [theta,rho] = cart2pol(intra_disc_loops{i}{j,1},intra_disc_loops{i}{j,2});
            rho=rho+abs(2*(min_axis_chr-sqrt((intra_disc_loops{i}{j,1}).^2 + (intra_disc_loops{i}{j,2}).^2)));
            [intra_disc_loops{i}{j,1},intra_disc_loops{i}{j,2}] = pol2cart(theta,rho);
        end
    end
    %resolved problem where two intertior parts of the loop are being reflected to opposite ends of the cylinder, the loop then crosses through the cylinder
    %reject walks with any steps of length spacing*2 or greater
    j=1;
    while j<=size(intra_disc_loops{i},1)-1
        while sqrt((intra_disc_loops{i}{j,1}-intra_disc_loops{i}{j+1,1}).^2 + (intra_disc_loops{i}{j,2}-intra_disc_loops{i}{j+1,2}).^2 + (intra_disc_loops{i}{j,3}-intra_disc_loops{i}{j+1,3}).^2) > spacing*2
            intra_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(intra_disc_loop_starts(i,:),intra_disc_loop_ends(i,:),round(intra_disc_loop_lengths(i)),spacing));
            j=1;
            for k=1:1:round(intra_disc_loop_lengths(i))
                if sqrt((intra_disc_loops{i}{k,1}).^2 + (intra_disc_loops{i}{k,2}).^2) <= min_axis_chr
                    [theta,rho] = cart2pol(intra_disc_loops{i}{k,1},intra_disc_loops{i}{k,2});
                    rho=rho+abs(2*(min_axis_chr-sqrt((intra_disc_loops{i}{k,1}).^2 + (intra_disc_loops{i}{k,2}).^2)));
                    [intra_disc_loops{i}{k,1},intra_disc_loops{i}{k,2}] = pol2cart(theta,rho);
                end
            end
        end
        j=j+1;
    end
    i
end

total_length_intra_disc = 0;
for i=1:1:num_intra_disc_loops
    total_length_intra_disc = total_length_intra_disc + round(intra_disc_loop_lengths(i));
end

% insert intra-disc loops into primary sequence in model
%%%RENAME NEW MODEL AND NEW_NEW_MODEL
new_new_model = insert_loops(num_intra_disc_loops, total_length_intra_disc, new_model, intra_disc_loops, intra_disc_loop_starts_ind);

figure
hold on
plot3(model(:,1),model(:,2),model(:,3))
for j = 1:1:size(inter_disc_loops,2)
    plot3(cell2mat(inter_disc_loops{j}(:,1)),cell2mat(inter_disc_loops{j}(:,2)),cell2mat(inter_disc_loops{j}(:,3)),'-','Color','r')
end
for k = 1:1:size(intra_disc_loops,2)
    plot3(cell2mat(intra_disc_loops{k}(:,1)),cell2mat(intra_disc_loops{k}(:,2)),cell2mat(intra_disc_loops{k}(:,3)),'-','Color','g')
end
xlim([min(new_model(:,1))*2 max(new_model(:,1))*2])
ylim([min(new_model(:,2))*2 max(new_model(:,2))*2])
zlim([min(new_model(:,3))*2 max(new_model(:,3))*2])
xlabel('x','FontSize', 24)
ylabel('y','FontSize', 24)
zlabel('z','FontSize', 24)
view(30,30)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numPoints = size(new_new_model,1);

MyColor = linspace(1,numPoints,numPoints)';
% create a connectivity matrix
Faces = [1:(numPoints-1); 2:numPoints]';

skip=1;

f=figure
hold on
plot3(new_new_model(:,1),new_new_model(:,2),new_new_model(:,3),'Color', [.6 .6 .6])
colormap jet
axis equal
xlim([min(model(:,1))*2 max(model(:,1))*2])
ylim([min(model(:,2))*2 max(model(:,2))*2])
zlim([min(model(:,3))*2 max(model(:,3))*2])
caxis([min(MyColor) max(MyColor)])
c = colorbar;
c.Label.String = 'primary sequence';
for i=skip+1:skip:numPoints
    clf(f)
    plot3(new_new_model(:,1),new_new_model(:,2),new_new_model(:,3),'Color', [.6 .6 .6])
    xlim([min(model(:,1))*2 max(model(:,1))*2])
    ylim([min(model(:,2))*2 max(model(:,2))*2])
    zlim([min(model(:,3))*2 max(model(:,3))*2])
    view(30,30)
    caxis([min(MyColor) max(MyColor)])
    c = colorbar;
    c.Label.String = 'primary sequence';
    patch('Faces', Faces(i-skip:i-1,:) ,'Vertices', new_new_model(1:i,:) ,'FaceColor', 'none', 'FaceVertexCData', MyColor(1:i,:) ,'EdgeColor','interp' ,'LineWidth',5, 'FaceAlpha',.5,'EdgeAlpha',.5);
    pause(0.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%y=y_0 + 1/(a sin(phi)*cos(phi))*log(abs(sin(x*a*sin(phi))+btan(phi))
%%%

%euclidian distance between monomers
%D = pdist(model); %in microns
D = pdist(new_model); %in microns
%D = pdist(new_new_model); %in microns

D = squareform(D);

%HiC contact probablities from distance
%P = spacing./D;
P = 1./D.^4;
P(isinf(P)) = 1; %self contact probabilities are 1

%are these coming from disc/loop connections?
%this is very important! as long extrachromosomal loops might have
%overlapping monomers
P(P>1) = 1; %no contact probabilities above 1

%output for GEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%writematrix(round(P,3,"significant"),'cholesteric_GEM_HiC.txt','Delimiter','tab')

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
        %PCSynth(count,1)=i*5000;
        %PCSynth(count,2)=j*5000;
        PCSynth(count,3)=P(i,j);
    end
end

PCSynth(:,3) = round(PCSynth(:,3),3,"significant");
writematrix(PCSynth,'cholesteric_CSynth_D4.txt','Delimiter','tab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot hic matrix
figure
imagesc(P);
colorbar
xlabel('~5kbp monomer index', 'fontsize', 24) %fix this!!!
ylabel('~5kbp monomer index', 'fontsize', 24)
%fix ticks!!!
set(gca,'ColorScale','log')
maj_axis_chr=colorbar;
maj_axis_chr.Label.String = 'Monomer Contact Probability';
maj_axis_chr.FontSize = 18;
%what are the boxes in this plot???

%total number of monomers
tot_monomers_in_body = sum(mon_per_layer);

%(The possible number of pairs of genomic positions separated by d on a given chromosome is Lc-d, where Lc is the length of the chromosome.)
distance_threshold = spacing*5;
%[s Ps] = contact_probability_xyz(model,distance_threshold);
% [s Ps] = contact_probability_xyz(new_model,distance_threshold);
%
% figure
% hold on
% ax=gca;
% xlim([0,400]);
% plot(s,Ps)
% legendStrings = "d_{cutoff} = " + string(distance_threshold) + "um";
% legend(legendStrings, 'Location', 'southwest')
% xlabel('s in multiples of 5kbp')
% ylabel('P(s)')
% set(ax,'xScale', 'log')
% set(ax,'YScale', 'log')
% grid on

function sum = nat_sum(x)
sum=0;
for i=1:x
    sum=sum+i;
end
end


function model_with_loops = insert_loops(numloops, total_length, model, loops, where_loops_start)
% Preallocate output
model_with_loops = zeros(size(model,1)+total_length,3);
% Find indices for old data
addRows = double(ismember(1:size(model,1), where_loops_start));

temp=find(addRows>0);
for i=1:1:numloops
    addRows(temp(i)) = size(loops{i},1);
end

oldDataInd = (1:size(model,1)) + cumsum([0, addRows(1:end-1)]);

% Add in old data
model_with_loops(oldDataInd,:) = model(:,:);

loops_concatenated=[];
for i=1:1:numloops
    loops_concatenated = [loops_concatenated; cell2mat(loops{i})];
end
% Find indices for new data
% Add in new data
model_with_loops(find(all(model_with_loops == 0,2)),:) = loops_concatenated;
end
