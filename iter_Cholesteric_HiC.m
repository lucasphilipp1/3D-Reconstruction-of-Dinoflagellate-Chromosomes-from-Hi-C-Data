clc
clear
%Author: Lucas Philipp
%REQUIRES STATISTICS AND MACHINE LEARNING TOOLBOX!

%%%% TO DO:
%%%% does a CLC model with a particular parameter set give rise to a unique scaling exponent? or can one scaling exponent correspond to multiple geometries?

%%% compute P(s) curve directly from HiC map, using cooltools without distance cutoff

num_chroms = 1;
s_avg = cell(1, num_chroms); %cell array for separation
Ps_avg = cell(1, num_chroms); %cell array for contact probability

for w = 1:1:num_chroms

    chromosome=[]; %position of DNA fibres
    %7.655Mbp/(10 layers)∗1monomer/5kbp≈(150 monomers)/layer
    num_mon = 200; %number of monomers in the thickest cholesteric disk

    frac_tot_sequence_in_loops = 0.4; %Fraction
    frac_loop_sequence_inter_disc = 0.5; %Fraction
    total_chromosome_length = 3000; %number of monomers

    %cholesteric pitch P: length in microns along long axis of the chromosome
    %corresponding to a full turn of the nematic director
    pitch = 1.30/2; %in microns

    %ellipse parameters for overall chromosome profile
    min_axis_chr=0.5; %radius of cholesteric disc, in microns

    %ellipse minor/major axes ratio for a cholesteric disc, for circle set ell_ratio=1
    %major axis of cholesteric disc ellipse = minor axis of chromosome ellipse
    ell_ratio = 0.8;

    spacing = min_axis_chr*sqrt(pi/num_mon); %mesh spacing in microns
    %2D grid of points
    x = -1:spacing:1; %mesh in microns
    y = -1:spacing:1; %mesh in microns

    %VECTOR FIELD ONLY:
    %longitudinal component magnitude
    b = 0.130 * 0.025; %in microns

    dtheta_layer=360*spacing/pitch; %DNA fibres in a layer are rotated dtheta_layer degrees counterclockwise relative to the previous layer

    % Find the target total length of core and loops
    wanted_total_core_length = total_chromosome_length-(frac_tot_sequence_in_loops*total_chromosome_length);
    wanted_total_inter_disc_loop_length = (total_chromosome_length-wanted_total_core_length)*frac_loop_sequence_inter_disc;
    wanted_total_intra_disc_loop_length = total_chromosome_length-wanted_total_core_length-wanted_total_inter_disc_loop_length;

    chol_layers = round(wanted_total_core_length/num_mon);
    mon_per_disc=zeros(chol_layers,1);
    maj_axis_chr=chol_layers*spacing; %in microns

    while size(chromosome,1) < wanted_total_core_length

        chromosome = [];
        discs = {};

        chol_layers = chol_layers +1;

        %for vector field describing orientation of DNA fibres
        lin_vec_x =[];
        lin_vec_y =[];

        intra_disc_loop_potential_starts_ind = [];

        i=1;
        for z = 0-(spacing*(chol_layers-1))/2:spacing:spacing*(chol_layers-1)-(spacing*(chol_layers-1))/2 %thickest/middle part of chromosome is z=0

            [X,Y] = meshgrid(x,y);

            % MAKE r A CONSTANT = maj_axis_chr (get longest distance parabola from origin)
            r = (min_axis_chr/maj_axis_chr)*sqrt(maj_axis_chr^2-z^2); %major axis of disk at this z position

            %mask for a single cholesteric disc
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
            intra_disc_loop_potential_starts_ind = [intra_disc_loop_potential_starts_ind; size(chromosome,1) + flip_these];

            chromosome = [chromosome; [rotcoord(:,1), rotcoord(:,2), z.*ones(size(xEllipse))]]; %append positions of DNA fibres for new cholesteric disc

            discs{end+1} = [rotcoord(:,1), rotcoord(:,2), z.*ones(size(xEllipse))];

            lin_vec_x = [lin_vec_x; spacing*cosd((i-1)*dtheta_layer).*ones(size(xEllipse))]; %append orientation of DNA fibres for new cholesteric disc
            lin_vec_y = [lin_vec_y; spacing*sind((i-1)*dtheta_layer).*ones(size(xEllipse))]; %append orientation of DNA fibres for new cholesteric disc

            mon_per_disc(i)=size(xEllipse,1);
            i=i+1;
        end

    end

    %orientation of DNA fibres, no longituindal component: fvec=(cos(az), sin(az), 0)
    lin_vec_0=ones(size(chromosome));
    %orientation of DNA fibres, longitudnal component=b fvec=(cos(az), sin(az), b) %normalize to unit length???
    lin_vec_b=ones(size(chromosome));

    lin_vec_0(:,1) = lin_vec_x;
    lin_vec_0(:,2) = lin_vec_y;
    lin_vec_0(:,3) = zeros(size(chromosome(:,3)));

    lin_vec_b(:,1) = lin_vec_x;
    lin_vec_b(:,2) = lin_vec_y;
    lin_vec_b(:,3) = b.*ones(size(chromosome(:,3))); %not affected by normalization

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
    chromosome_rot=chromosome*Rphi;
    lin_vec_0_rot=lin_vec_0*Rphi;
    lin_vec_b_rot=lin_vec_b*Rphi;
    % figure
    % set(gcf, 'Position', get(0, 'Screensize'));
    % q0=quiver3(chromosome_rot(:,1)-lin_vec_0_rot(:,1)./2,chromosome_rot(:,2)-lin_vec_0_rot(:,2)./2,chromosome_rot(:,3)-lin_vec_0_rot(:,3)./2,lin_vec_0_rot(:,1),lin_vec_0_rot(:,2),lin_vec_0_rot(:,3),'off'); %orientation vector rotates about midpoint instead of base
    % q0.ShowArrowHead = 'off';
    % xlim([-1 1])
    % ylim([-1 1])
    % zlim([-0.15 0.15])
    % view(0,90)
    % xlabel('x','FontSize', 24)
    % ylabel('y','FontSize', 24)
    % zlabel('z','FontSize', 24)
    % title('Rotated model for oblique cross-section viewing','FontSize', 24)

    %INTRA DISK LOOPS
    %choose potential intra-disc loop start positions, i.e. perimeter points on discs
    intra_disc_loop_potential_starts_ind = sort(intra_disc_loop_potential_starts_ind);

    append=[];
    %RENAME IDX2
    for i = 1:1:size(intra_disc_loop_potential_starts_ind,1)
        [match, idx2] = ismember(chromosome(intra_disc_loop_potential_starts_ind(i),:), chromosome, 'rows');
        intra_disc_loop_potential_starts_ind(i)=idx2;
        append = [append; idx2+1];
        append = [append; idx2-1];
    end
    intra_disc_loop_potential_starts_ind = [intra_disc_loop_potential_starts_ind; append];

    intra_disc_loop_potential_starts_ind = unique(intra_disc_loop_potential_starts_ind,'rows');

    intra_disc_loop_potential_starts_ind = sort(intra_disc_loop_potential_starts_ind);

    intra_disc_loop_potential_starts_ind(intra_disc_loop_potential_starts_ind == 0, :) = [];
    intra_disc_loop_potential_starts_ind(intra_disc_loop_potential_starts_ind == size(chromosome,1)+1, :) = [];

    save=[];
    for i = 3:1:size(intra_disc_loop_potential_starts_ind,1)-2
        if ismember(chromosome(intra_disc_loop_potential_starts_ind(i)-2,:), chromosome(intra_disc_loop_potential_starts_ind,:), 'rows') & ismember(chromosome(intra_disc_loop_potential_starts_ind(i)+2,:), chromosome(intra_disc_loop_potential_starts_ind,:), 'rows')
            save=[save; i];
        end
    end

    intra_disc_loop_potential_starts_ind(save)=[];

    % figure
    % hold on
    % plot3(chromosome(:,1),chromosome(:,2),chromosome(:,3),'-','Color','k')
    % plot3(chromosome(intra_disc_loop_potential_starts_ind,1),chromosome(intra_disc_loop_potential_starts_ind,2),chromosome(intra_disc_loop_potential_starts_ind,3),'o','MarkerFaceColor','r')
    % xlim([min(chromosome(:,1))*2 max(chromosome(:,1))*2])
    % ylim([min(chromosome(:,2))*2 max(chromosome(:,2))*2])
    % for i = 0:0.05:1
    % zlim([max(chromosome(:,3))*(i-0.04) max(chromosome(:,3))*(i+0.04)])
    % pause(5)
    % end
    % view(0,90)
    % hold off

    %INTER DISK LOOPS
    idx=find(diff(chromosome(:,3))>0); %new disc starts when z coordinate changes
    num_inter_disc_loops = size(idx,1);
    mean_looplength = wanted_total_inter_disc_loop_length/num_inter_disc_loops;
    % Generate inter disc loops with variable looplength drawn from an exponential distribution
    % Total length of inter disc loops adds up to the wanted length of inter disc loops
    inter_disc_loop_lengths = [];

    for i = 1:1:num_inter_disc_loops
        new_loop = round(exprnd(mean_looplength));
        % Reject loops lengths of 0
        while new_loop == 0
            new_loop = round(exprnd(mean_looplength));
        end
        inter_disc_loop_lengths = [inter_disc_loop_lengths new_loop];
    end

    inter_disc_loops = cell(1, num_inter_disc_loops); %cell array containing the loops

    disc_spine_1 = zeros(num_inter_disc_loops+1, 3); %on lattice points
    disc_spine_2 = zeros(num_inter_disc_loops+1, 3);  %on lattice points
    save_disc_starts = zeros(num_inter_disc_loops+1, 3); %remains off-latice, matches chromosome coordinates
    save_disc_ends = zeros(num_inter_disc_loops+1, 3); %remains off-latice, matches chromosome coordinates

    %find starts and ends to each disc
    for i=1:1:size(idx,1)+1
        if i == 1
            disc_spine_1(i,:)=chromosome(1,:);
            disc_spine_2(i,:)=chromosome(idx(1),:);
        elseif i == size(idx,1)+1
            disc_spine_1(i,:)=chromosome(idx(end)+1,:);
            disc_spine_2(i,:)=chromosome(end,:);
        else
            disc_spine_1(i,:)=chromosome(idx(i-1)+1,:);
            disc_spine_2(i,:)=chromosome(idx(i),:);
        end
    end

    for i=1:1:num_inter_disc_loops+1
        save_disc_starts(i,:)=disc_spine_1(i,:);
        save_disc_ends(i,:)=disc_spine_2(i,:);

        disc_spine_1(i,:)=spacing*round(disc_spine_1(i,:)./spacing);
        disc_spine_2(i,:)=spacing*round(disc_spine_2(i,:)./spacing);
    end

    count = 1;
    for i=1:1:num_inter_disc_loops
        i
        inter_disc_loop_lengths(i)
        %adjacent disks connected
        if mod(count,2)==1
            %convert 2D array to cell
            inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_spine_1(count,:),disc_spine_1(count+1,:),inter_disc_loop_lengths(i),spacing));
        else
            %convert 2D array to cell
            inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_spine_2(count,:),disc_spine_2(count+1,:),inter_disc_loop_lengths(i),spacing));
        end
        for j=1:1:inter_disc_loop_lengths(i)
            if sqrt((inter_disc_loops{i}{j,1}).^2 + (inter_disc_loops{i}{j,2}).^2) <= min_axis_chr*0.9
                [theta,rho] = cart2pol(inter_disc_loops{i}{j,1},inter_disc_loops{i}{j,2});
                rho=rho+abs(2*(min_axis_chr-sqrt((inter_disc_loops{i}{j,1}).^2 + (inter_disc_loops{i}{j,2}).^2)));
                [inter_disc_loops{i}{j,1},inter_disc_loops{i}{j,2}] = pol2cart(theta,rho);
            end
        end
        %resolved problem where two intertior parts of the loop are being reflected to opposite ends of the cylinder, the loop then crosses through the cylinder
        %reject walks with any steps of length spacing*3 or greater
        j=1;
        while j<=size(inter_disc_loops{i},1)-1
            while sqrt((inter_disc_loops{i}{j,1}-inter_disc_loops{i}{j+1,1}).^2 + (inter_disc_loops{i}{j,2}-inter_disc_loops{i}{j+1,2}).^2 + (inter_disc_loops{i}{j,3}-inter_disc_loops{i}{j+1,3}).^2) > spacing*4
                if mod(count,2)==1
                    inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_spine_1(count,:),disc_spine_1(count+1,:),inter_disc_loop_lengths(i),spacing));
                    j=1;
                else
                    inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_spine_2(count,:),disc_spine_2(count+1,:),inter_disc_loop_lengths(i),spacing));
                    j=1;
                end
                for k=1:1:inter_disc_loop_lengths(i)
                    if sqrt((inter_disc_loops{i}{k,1}).^2 + (inter_disc_loops{i}{k,2}).^2) <= min_axis_chr*0.9
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    chromosome_unflipped = chromosome;

    %by flipping every second disc the start of the next disc will be directly above the end of the disc below
    chromosome(1:idx(1),:) = flipud(chromosome(1:idx(1),:));
    for i=1:1:size(idx,1)
        if mod(i,2)==0
            if i == size(idx,1)
                chromosome(idx(i)+1:end,:) = flipud(chromosome(idx(i)+1:end,:));
            else
                chromosome(idx(i)+1:idx(i+1),:) = flipud(chromosome(idx(i)+1:idx(i+1),:));
            end
        end
    end

    for i = 1:1:size(intra_disc_loop_potential_starts_ind,1)
        [match, idx3] = ismember(chromosome_unflipped(intra_disc_loop_potential_starts_ind(i),:), chromosome, 'rows');
        intra_disc_loop_potential_starts_ind(i)=idx3;
    end

    %delete inter disc loops until wanted_total_inter_disc_loop_length is reached
    %discs are otherwise connected with a direct vertical step
    sum_inter_disc_loops=0;
    for i = 1:1:size(inter_disc_loops,2)
        sum_inter_disc_loops = sum_inter_disc_loops+size(inter_disc_loops{i},1);
    end

    % Generate loop length while the desired total length is not hit
    while sum_inter_disc_loops >= wanted_total_inter_disc_loop_length
        temp=randperm(size(inter_disc_loops,2),1);
        inter_disc_loops(:,temp) = [];
        idx(temp)=[];
        if size(temp,1) > 0
            num_inter_disc_loops = num_inter_disc_loops - size(temp,1);
        end
        sum_inter_disc_loops = 0;
        for i = 1:1:size(inter_disc_loops,2)
            sum_inter_disc_loops = sum_inter_disc_loops+size(inter_disc_loops{i},1);
        end
    end

    sum_inter_disc_loops = 0;
    for i = 1:1:size(inter_disc_loops,2)
        sum_inter_disc_loops = sum_inter_disc_loops+size(inter_disc_loops{i},1);
    end

    % insert inter-disc loops into primary sequence in model
    chromosome_w_inter_disc_loops = insert_loops(num_inter_disc_loops, sum_inter_disc_loops, chromosome, inter_disc_loops, idx);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Xcyl,Ycyl,Zcyl] = cylinder(min_axis_chr);
    Zcyl = Zcyl*(max(chromosome(:,3))-min(chromosome(:,3)));
    Zcyl = Zcyl - max(chromosome(:,3));

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

    % figure
    % plot3(new_model(:,1),new_model(:,2),new_model(:,3))
    % xlim([min(model(:,1))*2 max(model(:,1))*2])
    % ylim([min(model(:,2))*2 max(model(:,2))*2])
    % zlim([min(model(:,3))*2 max(model(:,3))*2])
    % xlabel('x','FontSize', 24)
    % ylabel('y','FontSize', 24)
    % zlabel('z','FontSize', 24)
    % view(30,30)

    %INTRA DISK LOOPS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %RENAME IDX3
    %this maps the list of indicies for intra-disc loop starts from chromosome -> chromosome_w_inter_disc_loops
    for i = 1:1:size(intra_disc_loop_potential_starts_ind,1)
        [match, idx3] = ismember(chromosome(intra_disc_loop_potential_starts_ind(i),:), chromosome_w_inter_disc_loops, 'rows');
        intra_disc_loop_potential_starts_ind(i)=idx3;
    end

    % figure
    % hold on
    % plot3(chromosome_w_inter_disc_loops(:,1),chromosome_w_inter_disc_loops(:,2),chromosome_w_inter_disc_loops(:,3),'-','Color','k')
    % plot3(chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,1),chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,2),chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,3),'o','MarkerFaceColor','r')
    % xlim([min(chromosome_w_inter_disc_loops(:,1))*2 max(chromosome_w_inter_disc_loops(:,1))*2])
    % ylim([min(chromosome_w_inter_disc_loops(:,2))*2 max(chromosome_w_inter_disc_loops(:,2))*2])
    % for i = 0:0.05:1
    % zlim([max(chromosome_w_inter_disc_loops(:,3))*(i-0.04) max(chromosome_w_inter_disc_loops(:,3))*(i+0.04)])
    % pause(5)
    % end
    % view(0,90)
    % hold off

    %ensure intra-disc loops don't start/end at the same position as an inter-disc loop
    for i = 1:1:size(save_disc_starts,1)
        common = ismember(save_disc_starts(i,:),chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,:),'rows');
        if common == 1
            RowIdx = find(ismember(chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,:), save_disc_starts(i,:),'rows'));
            intra_disc_loop_potential_starts_ind(RowIdx)=[];
        end
        common = ismember(save_disc_ends(i,:),chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,:),'rows');
        if common == 1
            RowIdx = find(ismember(chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,:), save_disc_ends(i,:),'rows'));
            intra_disc_loop_potential_starts_ind(RowIdx)=[];
        end
    end

    % Generate loops with looplength as mean from an exponential distribution
    % Total length of loops adds up to the wanted length of intra disc loops
    intra_disc_loop_lengths = [];

    while sum(intra_disc_loop_lengths) <= wanted_total_intra_disc_loop_length
        new_loop = round(exprnd(mean_looplength));
        % Reject loops lengths of 0
        while new_loop == 0
            new_loop = round(exprnd(mean_looplength));
        end
        if sum(intra_disc_loop_lengths)+new_loop > wanted_total_intra_disc_loop_length
            break; %CHECK THIS!!!
        end
        intra_disc_loop_lengths = [intra_disc_loop_lengths new_loop];
    end

    if sum(intra_disc_loop_lengths) < wanted_total_intra_disc_loop_length
        intra_disc_loop_lengths = [intra_disc_loop_lengths wanted_total_intra_disc_loop_length-sum(intra_disc_loop_lengths)];
    end

    num_intra_disc_loops = size(intra_disc_loop_lengths,2);
    intra_disc_loops = cell(1, num_intra_disc_loops); %cell array

    %intra_disc_loop_potential_starts are indicies for chromosome_w_inter_disc_loops
    cluster=clusterdata(intra_disc_loop_potential_starts_ind,1);
    for i = 1:1:max(cluster)
        del = [];
        if size(find(cluster==i),1) > 2
            del = [del; i];
        end
        intra_disc_loop_potential_starts_ind(find(del)) = [];
    end

    % figure
    % hold on
    % plot3(chromosome_w_inter_disc_loops(:,1),chromosome_w_inter_disc_loops(:,2),chromosome_w_inter_disc_loops(:,3),'-','Color','k')
    % plot3(chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,1),chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,2),chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,3),'o','MarkerFaceColor','r')
    % xlim([min(chromosome_w_inter_disc_loops(:,1))*2 max(chromosome_w_inter_disc_loops(:,1))*2])
    % ylim([min(chromosome_w_inter_disc_loops(:,2))*2 max(chromosome_w_inter_disc_loops(:,2))*2])
    % for i = 0:0.05:1
    % zlim([max(chromosome_w_inter_disc_loops(:,3))*(i-0.04) max(chromosome_w_inter_disc_loops(:,3))*(i+0.04)])
    % pause(5)
    % end
    % view(0,90)
    % hold off

    %number of random intra disc loops starts and ends correspond to num_intra_disc_loops and are chosen randomly from all possible locations
    intra_disc_loop_potential_starts_ind=intra_disc_loop_potential_starts_ind(find(diff(intra_disc_loop_potential_starts_ind)==1));
    p = randperm(size(intra_disc_loop_potential_starts_ind,1),num_intra_disc_loops);
    intra_disc_loop_starts_ind = sort(intra_disc_loop_potential_starts_ind(p));
    intra_disc_loop_ends_ind = intra_disc_loop_starts_ind + 1;

    intra_disc_loop_starts = zeros(num_intra_disc_loops+1, 3); %on lattice points
    intra_disc_loop_ends = zeros(num_intra_disc_loops+1, 3);  %on lattice points

    for i=1:1:num_intra_disc_loops
        intra_disc_loop_starts(i,:)=spacing*round(chromosome_w_inter_disc_loops(intra_disc_loop_starts_ind(i),:)./spacing);
        intra_disc_loop_ends(i,:)=spacing*round(chromosome_w_inter_disc_loops(intra_disc_loop_ends_ind(i),:)./spacing);
    end

    % figure
    % hold on
    % plot3(chromosome_w_inter_disc_loops(:,1),chromosome_w_inter_disc_loops(:,2),chromosome_w_inter_disc_loops(:,3))
    % plot3(chromosome_w_inter_disc_loops(intra_disc_loop_starts_ind,1),chromosome_w_inter_disc_loops(intra_disc_loop_starts_ind,2),chromosome_w_inter_disc_loops(intra_disc_loop_starts_ind,3),'o','MarkerFaceColor','r')
    % plot3(chromosome_w_inter_disc_loops(intra_disc_loop_ends_ind,1),chromosome_w_inter_disc_loops(intra_disc_loop_ends_ind,2),chromosome_w_inter_disc_loops(intra_disc_loop_ends_ind,3),'o','MarkerFaceColor','g')
    % xlim([min(chromosome_w_inter_disc_loops(:,1))*2 max(chromosome_w_inter_disc_loops(:,1))*2])
    % ylim([min(chromosome_w_inter_disc_loops(:,2))*2 max(chromosome_w_inter_disc_loops(:,2))*2])
    % zlim([min(chromosome_w_inter_disc_loops(:,3))*2 max(chromosome_w_inter_disc_loops(:,3))*2])
    % xlabel('x','FontSize', 24)
    % ylabel('y','FontSize', 24)
    % zlabel('z','FontSize', 24)
    % view(30,30)

    %Generate intra-disc loops
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:1:size(intra_disc_loop_lengths,2)
        intra_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(intra_disc_loop_starts(i,:),intra_disc_loop_ends(i,:),intra_disc_loop_lengths(i),spacing));
        for j=1:1:intra_disc_loop_lengths(i)
            if sqrt((intra_disc_loops{i}{j,1}).^2 + (intra_disc_loops{i}{j,2}).^2) <= min_axis_chr*0.9
                [theta,rho] = cart2pol(intra_disc_loops{i}{j,1},intra_disc_loops{i}{j,2});
                rho=rho+abs(2*(min_axis_chr-sqrt((intra_disc_loops{i}{j,1}).^2 + (intra_disc_loops{i}{j,2}).^2)));
                [intra_disc_loops{i}{j,1},intra_disc_loops{i}{j,2}] = pol2cart(theta,rho);
            end
        end
        %resolved problem where two intertior parts of the loop are being reflected to opposite ends of the cylinder, the loop then crosses through the cylinder
        %reject walks with any steps of length spacing*3 or greater
        j=1;
        while j<=size(intra_disc_loops{i},1)-1
            while sqrt((intra_disc_loops{i}{j,1}-intra_disc_loops{i}{j+1,1}).^2 + (intra_disc_loops{i}{j,2}-intra_disc_loops{i}{j+1,2}).^2 + (intra_disc_loops{i}{j,3}-intra_disc_loops{i}{j+1,3}).^2) > spacing*4
                intra_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(intra_disc_loop_starts(i,:),intra_disc_loop_ends(i,:),intra_disc_loop_lengths(i),spacing));
                j=1;
                for k=1:1:intra_disc_loop_lengths(i)
                    if sqrt((intra_disc_loops{i}{k,1}).^2 + (intra_disc_loops{i}{k,2}).^2) <= min_axis_chr*0.9
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

    % insert intra-disc loops into primary sequence in chromosome
    chromosome_w_inter_and_intra_disc_loops = insert_loops(num_intra_disc_loops, sum(intra_disc_loop_lengths), chromosome_w_inter_disc_loops, intra_disc_loops, intra_disc_loop_starts_ind);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    numPoints = size(chromosome_w_inter_and_intra_disc_loops,1);

    MyColor = linspace(1,numPoints,numPoints)';
    % create a connectivity matrix
    Faces = [1:(numPoints-1); 2:numPoints]';

    skip=numPoints;

    f=figure
    hold on
    plot3(chromosome_w_inter_and_intra_disc_loops(:,1),chromosome_w_inter_and_intra_disc_loops(:,2),chromosome_w_inter_and_intra_disc_loops(:,3),'Color', [.6 .6 .6])
    colormap jet
    axis equal
    xlim([min(chromosome(:,1))*2 max(chromosome(:,1))*2])
    ylim([min(chromosome(:,2))*2 max(chromosome(:,2))*2])
    zlim([min(chromosome(:,3))*2 max(chromosome(:,3))*2])
    caxis([min(MyColor) max(MyColor)])
    c = colorbar;
    c.Label.String = 'primary sequence';
    patch('Faces', Faces(:,:) ,'Vertices', chromosome_w_inter_and_intra_disc_loops(:,:) ,'FaceColor', 'none', 'FaceVertexCData', MyColor(:,:) ,'EdgeColor','interp' ,'LineWidth',5, 'FaceAlpha',.5,'EdgeAlpha',.5);
    view(30,30)
    % for i=skip+1:skip:numPoints
    %     clf(f)
    %     plot3(new_new_model(:,1),new_new_model(:,2),new_new_model(:,3),'Color', [.6 .6 .6])
    %     xlim([min(model(:,1))*2 max(model(:,1))*2])
    %     ylim([min(model(:,2))*2 max(model(:,2))*2])
    %     zlim([min(model(:,3))*2 max(model(:,3))*2])
    %     view(30,30)
    %     caxis([min(MyColor) max(MyColor)])
    %     c = colorbar;
    %     c.Label.String = 'primary sequence';
    %     patch('Faces', Faces(i-skip:i-1,:) ,'Vertices', new_new_model(1:i,:) ,'FaceColor', 'none', 'FaceVertexCData', MyColor(1:i,:) ,'EdgeColor','interp' ,'LineWidth',5, 'FaceAlpha',.5,'EdgeAlpha',.5);
    %     pause(0.5);
    % end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%
    %y=y_0 + 1/(a sin(phi)*cos(phi))*log(abs(sin(x*a*sin(phi))+btan(phi))
    %%%

    %euclidian distance between monomers
    D = pdist(chromosome_w_inter_and_intra_disc_loops); %in microns

    D = squareform(D);

    %HiC contact probablities from distance
    %P = spacing./D;
    P = spacing./D.^4;
    P(isinf(P)) = 1; %self contact probabilities are 1

    %are these coming from disc/loop connections?
    %this is very important! as long extrachromosomal loops might have
    %overlapping monomers
    P(P>1) = 1; %no contact probabilities above 1

    % figure
    % hold on
    % plot3(chromosome(:,1),chromosome(:,2),chromosome(:,3))
    % for j = 1:1:size(inter_disc_loops,2)
    %     plot3(cell2mat(inter_disc_loops{j}(:,1)),cell2mat(inter_disc_loops{j}(:,2)),cell2mat(inter_disc_loops{j}(:,3)),'-','Color','r')
    % end
    % for k = 1:1:size(intra_disc_loops,2)
    %     plot3(cell2mat(intra_disc_loops{k}(:,1)),cell2mat(intra_disc_loops{k}(:,2)),cell2mat(intra_disc_loops{k}(:,3)),'-','Color','g')
    % end
    % for i = 1:1:size(find(diag(D,1)>0.5),1)
    %     far=find(diag(D,1)>0.5);
    %     plot3(chromosome_w_inter_and_intra_disc_loops(far(i):far(i)+1,1),chromosome_w_inter_and_intra_disc_loops(far(i):far(i)+1,2),chromosome_w_inter_and_intra_disc_loops(far(i):far(i)+1,3),'-o','Color','k')
    % end
    % xlim([min(chromosome_w_inter_disc_loops(:,1))*2 max(chromosome_w_inter_disc_loops(:,1))*2])
    % ylim([min(chromosome_w_inter_disc_loops(:,2))*2 max(chromosome_w_inter_disc_loops(:,2))*2])
    % zlim([min(chromosome_w_inter_disc_loops(:,3))*2 max(chromosome_w_inter_disc_loops(:,3))*2])
    % xlabel('x','FontSize', 24)
    % ylabel('y','FontSize', 24)
    % zlabel('z','FontSize', 24)
    % view(30,30)

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
    %writematrix(PCSynth,'cholesteric_CSynth_D4.txt','Delimiter','tab')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    discs_indicies = [];
    discs_starts = [];
    discs_ends = [];
    intra_disc_loop_indicies = [];

    for i = 1:1:size(discs,2)
        discs_starts = [discs_starts; min(find(ismember(chromosome_w_inter_and_intra_disc_loops, discs{i},'rows')))];
        discs_ends = [discs_ends; max(find(ismember(chromosome_w_inter_and_intra_disc_loops, discs{i},'rows')))];
        discs_indicies = [discs_indicies; [min(find(ismember(chromosome_w_inter_and_intra_disc_loops, discs{i},'rows'))):max(find(ismember(chromosome_w_inter_and_intra_disc_loops, discs{i},'rows')))]']; %highlight includes intra-disc loops
        intra_disc_loop_indicies = [intra_disc_loop_indicies;     setxor([min(find(ismember(chromosome_w_inter_and_intra_disc_loops, discs{i},'rows'))):max(find(ismember(chromosome_w_inter_and_intra_disc_loops, discs{i},'rows')))]', find(ismember(chromosome_w_inter_and_intra_disc_loops, discs{i},'rows')))]; 
    end

    f1=figure
    %plot cholesteric disc track
    hold on
    histogram(discs_indicies, size(chromosome_w_inter_and_intra_disc_loops,1))
    for i = 1:1:size(discs_ends,1)
        xline(discs_starts(i),'w','LineWidth',3)
        xline(discs_ends(i),'w','LineWidth',3)
    end
    title('Cholesteric Disc Locations')
    box off
    axis off
    daspect([100 1 1]) %squish height of figure
    hold off
    ax1=gca;

    % f2=figure
    % %plot intra-disc loops track
    % histogram(intra_disc_loop_indicies, size(chromosome_w_inter_and_intra_disc_loops,1), 'EdgeColor', 'r')
    % title('Intra Disc Loop Locations')
    % box off
    % axis off
    % daspect([100 1 1]) %squish height of figure
    % ax2=gca;

    f2=figure
    %plot hic matrix
    imagesc(P);
    hold on;
    colorbar
    xlabel('monomer #', 'fontsize', 24)
    ylabel('monomer #', 'fontsize', 24)
    %fix ticks!!!
    set(gca,'ColorScale','log')
    maj_axis_chr=colorbar;
    maj_axis_chr.Label.String = 'Monomer Contact Probability';
    maj_axis_chr.FontSize = 18;
    ax2=gca;

    figure
    tcl=tiledlayout(8,1); % 8 by 1 tile layout
    ax1.Parent=tcl;
    ax1.Layout.Tile=1; % placing fig1 in the 1st tile
    ax2.Parent=tcl;
    ax2.Layout.Tile=2; % placing fig2 in the 2nd tile
    % ax3.Parent=tcl;
    % ax3.Layout.Tile=3; % placing fig3 in the 3nd tile
    ax2.Layout.TileSpan=[7 1];  %fig2 spans tiles 2 through 8
    linkaxes([ax1 ax2],'x')
    xlim([0 size(chromosome_w_inter_and_intra_disc_loops,1)])
    close(f1)
    close(f2)

    %(The possible number of pairs of genomic positions separated by d on a given chromosome is Lc-d, where Lc is the length of the chromosome.)
    distance_threshold = spacing*2;
    [s Ps] = contact_probability_xyz(chromosome_w_inter_and_intra_disc_loops,distance_threshold);

    s_avg{w} = s';
    Ps_avg{w} = Ps';

end

% Vertically concatenate, pad with NaNs
maxNumCol = max(cellfun(@(c) size(c,2), Ps_avg));  % max number of columns
aMat = cell2mat(cellfun(@(c){padarray(c,[0,maxNumCol-size(c,2)],NaN,'Post')}, Ps_avg)');
colMeans = mean(aMat,1,'omitnan')';

Ps_avg = cellfun(@transpose,Ps_avg,'UniformOutput',false);
s_avg = cellfun(@transpose,s_avg,'UniformOutput',false);

[max_size, max_index] = max(cellfun('size', s_avg, 1));
[min_size, min_index] = min(cellfun('size', s_avg, 1));

figure
hold on
ax=gca;
err = [];
for v = 1:1:size(s_avg{min_index},1)
    spread_Ps = [];
    for w = 1:1:num_chroms
        spread_Ps = [spread_Ps Ps_avg{w}(v)];
    end
    err = [err std(spread_Ps)];
end
errorbar(s_avg{min_index},colMeans(1:size(s_avg{min_index},1)),err)
%legendStrings = "d_{cutoff} = " + string(distance_threshold) + "um";
%legend(legendStrings, 'Location', 'southwest')
xlabel('s [# of monomers]')
ylabel('P(s)')
set(ax,'xScale', 'log')
set(ax,'YScale', 'log')
xlim([0,total_chromosome_length]);
grid on

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