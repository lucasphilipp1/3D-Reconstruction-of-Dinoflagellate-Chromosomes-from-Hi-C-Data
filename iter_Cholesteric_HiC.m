clc
clear
%Author: Lucas Philipp
%iteratively create CLC chromosomes with various parameters
%simulate CLC Hi-C contact maps/contact probability curves
%
%requires: statistics and machine learning toolbox

num_chroms = 100; %number of simulated chromosomes
resolution = 5000; %number of base pairs per monomer

%contact probability curve data
s_avg = cell(1, num_chroms); %cell array for genomic separation
Ps_avg = cell(1, num_chroms); %cell array for contact probability

total_chromosome_length = 4000; %number of monomers
P_agg = zeros(total_chromosome_length); %collect Hi-C contact maps for different chromosomes

disc_diameter = []; %collect disc_diameter for different chromosomes
mon_per_disc_agg = []; %collect number of monomers per disc for different chromosomes
mean_loop_length_agg = []; %collect mean loop length for different chromosomes

for d = 1:1:num_chroms
    d
    chromosome=[]; %position of DNA
    %7.655Mbp/(10 layers)∗1monomer/5kbp≈(150 monomers)/layer
    num_mon_per_disc = 200; %number of monomers in the thickest cholesteric disc
    %num_mon_per_disc = random('Normal',150,40); %number of monomers in the thickest cholesteric disc
    %num_mon_per_disc = random('Uniform',75,150); %number of monomers in the thickest cholesteric disc

    frac_tot_sequence_in_loops = 0.1; %Fraction
    frac_loop_sequence_inter_disc = 0.5; %Fraction

    %cholesteric pitch P: length in microns along long axis of the chromosome
    %corresponding to a full turn of the nematic director
    pitch = 1.30/2; %in microns

    %ellipse parameters for overall chromosome profile
    min_axis_chr=0.5; %radius of cholesteric disc, in microns

    % optional slanted discs (discs not perpendicular to chromosome long axis)
    %tilt_angle = 10; %in degrees, 0 is perpendicular to chromosome long axis
    %tilt_angle = random('Normal',0,30); %in degrees

    xy_spacing = min_axis_chr*sqrt(pi/num_mon_per_disc); %mesh spacing in microns
    z_spacing = xy_spacing;

    %distance cutoff to count as HiC contact
    distance_threshold = xy_spacing*2;

    %2D grid of points
    x = -1:xy_spacing:1; %mesh in microns
    y = -1:xy_spacing:1; %mesh in microns

    %VECTOR FIELD ONLY:
    %longitudinal component magnitude
    b = 0.130 * 0.025; %in microns

    dtheta_layer=360*z_spacing/pitch; %DNA fibres in a layer are rotated dtheta_layer degrees counterclockwise relative to the previous layer

    % Find the target total length of sequence in core and in loops
    wanted_total_core_length = total_chromosome_length-(frac_tot_sequence_in_loops*total_chromosome_length);
    wanted_total_inter_disc_loop_length = (total_chromosome_length-wanted_total_core_length)*frac_loop_sequence_inter_disc;
    wanted_total_intra_disc_loop_length = total_chromosome_length-wanted_total_core_length-wanted_total_inter_disc_loop_length;

    chol_layers = round(wanted_total_core_length/num_mon_per_disc);
    mon_per_disc=zeros(chol_layers,1);
    maj_axis_chr=chol_layers*z_spacing; %in microns

    while size(chromosome,1) < wanted_total_core_length

        chromosome = [];
        %discs are different sizes, need cell array
        discs = {};

        chol_layers = chol_layers +1;

        %for vector field describing orientation of DNA fibres
        lin_vec_x =[];
        lin_vec_y =[];

        intra_disc_loop_potential_starts_ind = [];

        i=1;
        for z = 0-(z_spacing*(chol_layers-1))/2:z_spacing:z_spacing*(chol_layers-1)-(z_spacing*(chol_layers-1))/2 %thickest/middle part of chromosome is z=0

            [X,Y] = meshgrid(x,y);

            %chromosome discs are wider at center than at either end
            r = (min_axis_chr/maj_axis_chr)*sqrt(maj_axis_chr^2-z^2); %major axis of disc at this z position

            %mask for a single cholesteric disc
            indicator = sqrt((X./r).^2 + (Y./r).^2) - 1 < 0;
            xEllipse = X(indicator);
            yEllipse = Y(indicator);

            %every second line of monomers y->-y so ordering does a "weave"
            [c,ia,ib] = unique(xEllipse);
            flip_these=[];
            max_diameter = 0;
            for k = 1:2:length(c)
                flip_these = [flip_these; find(ib==k)];
                temp = size(find(ib==k),1);
                if temp > max_diameter
                    max_diameter=temp;
                end
            end

            disc_diameter = [disc_diameter max_diameter];

            yEllipse(flip_these) = -1.*yEllipse(flip_these);

            %perform rotation for each disc
            R=[cosd((i-1)*dtheta_layer) -sind((i-1)*dtheta_layer); sind((i-1)*dtheta_layer) cosd((i-1)*dtheta_layer)]; %create 2D rotation matix, applies to (X,Y) pairs
            rotcoord=[];
            for j=1:1:size(xEllipse,1)
                rotcoord=[rotcoord; [xEllipse(j), yEllipse(j)]*R']; %multiply by rotation matrix
            end

            %create array of potential start positions for intra-disc extra chromosomal loops
            intra_disc_loop_potential_starts_ind = [intra_disc_loop_potential_starts_ind; size(chromosome,1) + flip_these];

            chromosome = [chromosome; [rotcoord(:,1), rotcoord(:,2), z.*ones(size(xEllipse))]]; %append positions of DNA fibres for new cholesteric disc

            discs{end+1} = [rotcoord(:,1), rotcoord(:,2), z.*ones(size(xEllipse))];

            lin_vec_x = [lin_vec_x; xy_spacing*cosd((i-1)*dtheta_layer).*ones(size(xEllipse))]; %append orientation of DNA fibres for new cholesteric disc
            lin_vec_y = [lin_vec_y; xy_spacing*sind((i-1)*dtheta_layer).*ones(size(xEllipse))]; %append orientation of DNA fibres for new cholesteric disc

            mon_per_disc(i)=size(xEllipse,1);
            i=i+1;
        end

    end
    mon_per_disc_agg = [mon_per_disc_agg; mon_per_disc];

    %%%
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

    %visualize oblique TEM cross-section through CLC model
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
    %%%

    %INTRA DISC LOOPS
    %choose potential intra-disc loop start positions, i.e. nucleoplasm-exposed monomers on discs
    intra_disc_loop_potential_starts_ind = sort(intra_disc_loop_potential_starts_ind);

    append=[];
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

    % visualize intra-disc loop potential start positions
    % figure
    % hold on
    % plot3(chromosome(:,1),chromosome(:,2),chromosome(:,3),'-','Color','k')
    % plot3(chromosome(intra_disc_loop_potential_starts_ind,1),chromosome(intra_disc_loop_potential_starts_ind,2),chromosome(intra_disc_loop_potential_starts_ind,3),'o','MarkerFaceColor','r')
    % xlim([min(chromosome(:,1))*2 max(chromosome(:,1))*2])
    % ylim([min(chromosome(:,2))*2 max(chromosome(:,2))*2])
    % set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
    % for i = 0:0.05:1
    %     zlim([max(chromosome(:,3))*(i-0.03) max(chromosome(:,3))*(i+0.03)])
    %     pause(5)
    % end
    % view(0,90)
    % hold off

    %INTER DISC LOOPS
    idx=find(diff(chromosome(:,3))>0); %new disc starts when z coordinate changes
    num_inter_disc_loops = size(idx,1);
    mean_loop_length = wanted_total_inter_disc_loop_length/num_inter_disc_loops;
    mean_loop_length_agg = [mean_loop_length_agg; mean_loop_length];
    % Generate inter disc loops with variable looplength drawn from an exponential distribution
    % Total length of inter disc loops adds up to the wanted length of inter disc loops
    inter_disc_loop_lengths = [];

    for i = 1:1:num_inter_disc_loops
        new_loop = round(exprnd(mean_loop_length));
        % Reject loops lengths of 0
        while new_loop == 0
            new_loop = round(exprnd(mean_loop_length));
        end
        inter_disc_loop_lengths = [inter_disc_loop_lengths new_loop];
    end

    inter_disc_loops = cell(1, num_inter_disc_loops); %cell array containing the loops

    %because of disc rotation all chromosome monomers are no longer on the same lattice
    %loop random walks are on the original xyz lattice
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

        disc_spine_1(i,:)=xy_spacing*round(disc_spine_1(i,:)./xy_spacing);
        disc_spine_2(i,:)=xy_spacing*round(disc_spine_2(i,:)./xy_spacing);
    end

    %extra-chromosomal loops avoid the condensed chromosome core
    count = 1;
    for i=1:1:num_inter_disc_loops
        %adjacent discs connected
        if mod(count,2)==1
            %convert 2D array to cell
            inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_spine_1(count,:),disc_spine_1(count+1,:),inter_disc_loop_lengths(i),z_spacing));
        else
            %convert 2D array to cell
            inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_spine_2(count,:),disc_spine_2(count+1,:),inter_disc_loop_lengths(i),z_spacing));
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
            while sqrt((inter_disc_loops{i}{j,1}-inter_disc_loops{i}{j+1,1}).^2 + (inter_disc_loops{i}{j,2}-inter_disc_loops{i}{j+1,2}).^2 + (inter_disc_loops{i}{j,3}-inter_disc_loops{i}{j+1,3}).^2) > z_spacing*4
                if mod(count,2)==1
                    inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_spine_1(count,:),disc_spine_1(count+1,:),inter_disc_loop_lengths(i),z_spacing));
                    j=1;
                else
                    inter_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(disc_spine_2(count,:),disc_spine_2(count+1,:),inter_disc_loop_lengths(i),z_spacing));
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

    %generate another loop while the desired total amount of sequence in loops is not yet reached
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

    % insert inter-disc loops into proper primary sequence locations in 3D model
    chromosome_w_inter_disc_loops = insert_loops(num_inter_disc_loops, sum_inter_disc_loops, chromosome, inter_disc_loops, idx);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %INTRA DISC LOOPS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %this maps the list of monomer positions for intra-disc loop starts from chromosome -> chromosome_w_inter_disc_loops
    for i = 1:1:size(intra_disc_loop_potential_starts_ind,1)
        [match, idx3] = ismember(chromosome(intra_disc_loop_potential_starts_ind(i),:), chromosome_w_inter_disc_loops, 'rows');
        intra_disc_loop_potential_starts_ind(i)=idx3;
    end

    % visualize that intra-disc loop coordinates were properly mapped
    % figure
    % hold on
    % plot3(chromosome_w_inter_disc_loops(:,1),chromosome_w_inter_disc_loops(:,2),chromosome_w_inter_disc_loops(:,3),'o','MarkerFaceColor','k')
    % plot3(chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,1),chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,2),chromosome_w_inter_disc_loops(intra_disc_loop_potential_starts_ind,3),'o','MarkerFaceColor','r')
    % xlim([min(chromosome_w_inter_disc_loops(:,1))*2 max(chromosome_w_inter_disc_loops(:,1))*2])
    % ylim([min(chromosome_w_inter_disc_loops(:,2))*2 max(chromosome_w_inter_disc_loops(:,2))*2])
    % set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
    % i=0;
    % for i = 0:0.05:1
    %     zlim([min(chromosome_w_inter_disc_loops(:,3))+(i-0.03) min(chromosome_w_inter_disc_loops(:,3))+(i+0.03)])
    %    pause(5)
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

    %generate loops with looplength as mean from an exponential distribution
    %total length of loops adds up to the wanted length of intra disc loops
    intra_disc_loop_lengths = [];

    while sum(intra_disc_loop_lengths) <= wanted_total_intra_disc_loop_length
        new_loop = round(exprnd(mean_loop_length));
        % Reject loops lengths of 0
        while new_loop == 0
            new_loop = round(exprnd(mean_loop_length));
        end
        if sum(intra_disc_loop_lengths)+new_loop > wanted_total_intra_disc_loop_length
            break;
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

    %number of random intra disc loops starts and ends correspond to num_intra_disc_loops and are chosen randomly from all possible locations
    intra_disc_loop_potential_starts_ind=intra_disc_loop_potential_starts_ind(find(diff(intra_disc_loop_potential_starts_ind)==1));
    p = randperm(size(intra_disc_loop_potential_starts_ind,1),num_intra_disc_loops);
    intra_disc_loop_starts_ind = sort(intra_disc_loop_potential_starts_ind(p));
    intra_disc_loop_ends_ind = intra_disc_loop_starts_ind + 1;

    intra_disc_loop_starts = zeros(num_intra_disc_loops+1, 3); %on lattice points
    intra_disc_loop_ends = zeros(num_intra_disc_loops+1, 3);  %on lattice points

    for i=1:1:num_intra_disc_loops
        intra_disc_loop_starts(i,:)=xy_spacing*round(chromosome_w_inter_disc_loops(intra_disc_loop_starts_ind(i),:)./xy_spacing);
        intra_disc_loop_ends(i,:)=xy_spacing*round(chromosome_w_inter_disc_loops(intra_disc_loop_ends_ind(i),:)./xy_spacing);
    end

    %visualize intra-disc loop start/end positions
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

    %generate intra-disc loops
    %random walk avoids condensed chromosome core
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:1:size(intra_disc_loop_lengths,2)
        intra_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(intra_disc_loop_starts(i,:),intra_disc_loop_ends(i,:),intra_disc_loop_lengths(i),z_spacing));
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
            while sqrt((intra_disc_loops{i}{j,1}-intra_disc_loops{i}{j+1,1}).^2 + (intra_disc_loops{i}{j,2}-intra_disc_loops{i}{j+1,2}).^2 + (intra_disc_loops{i}{j,3}-intra_disc_loops{i}{j+1,3}).^2) > z_spacing*4
                intra_disc_loops{i}=num2cell(constrained_self_avoiding_RW_3D(intra_disc_loop_starts(i,:),intra_disc_loop_ends(i,:),intra_disc_loop_lengths(i),z_spacing));
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
    end

    % insert intra-disc loops into primary sequence in chromosome
    chromosome_w_inter_and_intra_disc_loops = insert_loops(num_intra_disc_loops, sum(intra_disc_loop_lengths), chromosome_w_inter_disc_loops, intra_disc_loops, intra_disc_loop_starts_ind);

    % optional slanted discs (discs not perpendicular to chromosome long axis)
    tilt_bins = 50;
    if abs(tilt_angle) > 0
        tilt_x_coords=linspace(min(chromosome_w_inter_and_intra_disc_loops(:,1)),max(chromosome_w_inter_and_intra_disc_loops(:,1)),tilt_bins);
        tilt_z_coords=linspace(0,atan(tilt_angle*pi/180)*(max(tilt_x_coords)-min(tilt_x_coords)),tilt_bins-1);
        for i = 1:1:tilt_bins-1
            chromosome_w_inter_and_intra_disc_loops(find(tilt_x_coords(i+1) > chromosome_w_inter_and_intra_disc_loops(:,1) & chromosome_w_inter_and_intra_disc_loops(:,1) > tilt_x_coords(i)),3) = chromosome_w_inter_and_intra_disc_loops(find(tilt_x_coords(i+1) > chromosome_w_inter_and_intra_disc_loops(:,1) & chromosome_w_inter_and_intra_disc_loops(:,1) > tilt_x_coords(i)),3) + tilt_z_coords(i);
        end
    else
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % visualize primary sequence throughout chromosome
    numPoints = size(chromosome_w_inter_and_intra_disc_loops,1);

    %low res CLC structure
    %chromosome_w_inter_and_intra_disc_loops_100kb_res = chromosome_w_inter_and_intra_disc_loops(1:20:end,:);
    %numPoints = size(chromosome_w_inter_and_intra_disc_loops_100kb_res,1);

    MyColor = linspace(1,numPoints,numPoints)';
    % create a connectivity matrix
    Faces = [1:(numPoints-1); 2:numPoints]';

    skip=numPoints;

    % f=figure
    % screen = get(0, 'Screensize');
    % screen(3)=screen(3)/1.75;
    % set(gcf, 'Position', screen);
    % hold on
    % plot3(chromosome_w_inter_and_intra_disc_loops(:,1),chromosome_w_inter_and_intra_disc_loops(:,2),chromosome_w_inter_and_intra_disc_loops(:,3),'Color', [.6 .6 .6])
    % %plot3(chromosome_w_inter_and_intra_disc_loops_100kb_res(:,1),chromosome_w_inter_and_intra_disc_loops_100kb_res(:,2),chromosome_w_inter_and_intra_disc_loops_100kb_res(:,3),'Color', [.6 .6 .6]) %low res CLC structure
    % colormap jet
    % axis equal
    % xlim([min(chromosome_w_inter_and_intra_disc_loops(:,1))*1.5 max(chromosome_w_inter_and_intra_disc_loops(:,1))*1.5])
    % ylim([min(chromosome_w_inter_and_intra_disc_loops(:,2))*1.5 max(chromosome_w_inter_and_intra_disc_loops(:,2))*1.5])
    % zlim([min(chromosome_w_inter_and_intra_disc_loops(:,3))*1.1 max(chromosome_w_inter_and_intra_disc_loops(:,3))*1.1])
    % set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
    % caxis([min(MyColor) max(MyColor)])
    % c = colorbar;
    % c.Position = c.Position - [.1 0 0 0];
    % c.Ticks = linspace(0, total_chromosome_length, round(total_chromosome_length/500)+1);
    % c.TickLabels = num2cell(linspace(0, total_chromosome_length*resolution, round(total_chromosome_length/500)+1));
    % c.Label.String = 'primary sequence [bp]';
    % c.FontSize = 32;
    % patch('Faces', Faces(:,:) ,'Vertices', chromosome_w_inter_and_intra_disc_loops(:,:) ,'FaceColor', 'none', 'FaceVertexCData', MyColor(:,:) ,'EdgeColor','interp' ,'LineWidth',5, 'FaceAlpha',.5,'EdgeAlpha',.5);
    % %patch('Faces', Faces(:,:) ,'Vertices', chromosome_w_inter_and_intra_disc_loops_100kb_res(:,:) ,'FaceColor', 'none', 'FaceVertexCData', MyColor(:,:) ,'EdgeColor','interp' ,'LineWidth',5, 'FaceAlpha',.5,'EdgeAlpha',.5); %low res CLC structure
    % view(90,0)

    % for i=skip+1:skip:numPoints
    %     clf(f)
    %     plot3(chromosome_w_inter_and_intra_disc_loops(:,1),chromosome_w_inter_and_intra_disc_loops(:,2),chromosome_w_inter_and_intra_disc_loops(:,3),'Color', [.6 .6 .6])
    %     xlim([min(chromosome_w_inter_and_intra_disc_loops(:,1))*2 max(chromosome_w_inter_and_intra_disc_loops(:,1))*2])
    %     ylim([min(chromosome_w_inter_and_intra_disc_loops(:,2))*2 max(chromosome_w_inter_and_intra_disc_loops(:,2))*2])
    %     zlim([min(chromosome_w_inter_and_intra_disc_loops(:,3))*2 max(chromosome_w_inter_and_intra_disc_loops(:,3))*2])
    %     view(30,30)
    %     caxis([min(MyColor) max(MyColor)])
    %     c = colorbar;
    %     c.Label.String = 'primary sequence';
    %     patch('Faces', Faces(i-skip:i-1,:) ,'Vertices', chromosome_w_inter_and_intra_disc_loops(1:i,:) ,'FaceColor', 'none', 'FaceVertexCData', MyColor(1:i,:) ,'EdgeColor','interp' ,'LineWidth',5, 'FaceAlpha',.5,'EdgeAlpha',.5);
    %     pause(0.5);
    % end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %simluate Hi-C contact map
    %euclidian distance between monomers
    D = pdist(chromosome_w_inter_and_intra_disc_loops); %in microns

    D = squareform(D);

    %HiC contact probablities from distance
    %P = spacing./D;
    P = (xy_spacing/10)./D.^4;
    P(isinf(P)) = 1; %self contact probabilities are 1
    P(P>1) = 1; %no contact probabilities above 1

    %figure S1B in paper
    figure
    hold on
    plot3(chromosome(:,1),chromosome(:,2),chromosome(:,3))
    for j = 1:1:size(inter_disc_loops,2)
        plot3(cell2mat(inter_disc_loops{j}(:,1)),cell2mat(inter_disc_loops{j}(:,2)),cell2mat(inter_disc_loops{j}(:,3)),'-o','LineWidth',5,'Color','r')
    end
    for k = 1:1:size(intra_disc_loops,2)
        plot3(cell2mat(intra_disc_loops{k}(:,1)),cell2mat(intra_disc_loops{k}(:,2)),cell2mat(intra_disc_loops{k}(:,3)),'-o','LineWidth',5,'Color','g')
    end
    % for i = 1:1:size(find(diag(D,1)>0.5),1)
    %     far=find(diag(D,1)>0.5);
    %     plot3(chromosome_w_inter_and_intra_disc_loops(far(i):far(i)+1,1),chromosome_w_inter_and_intra_disc_loops(far(i):far(i)+1,2),chromosome_w_inter_and_intra_disc_loops(far(i):far(i)+1,3),'-o','Color','k')
    % end
    plot3(chromosome_w_inter_and_intra_disc_loops(:,1),chromosome_w_inter_and_intra_disc_loops(:,2),chromosome_w_inter_and_intra_disc_loops(:,3),'-','Color','k')
    xlim([min(chromosome_w_inter_disc_loops(:,1))*2 max(chromosome_w_inter_disc_loops(:,1))*2])
    ylim([min(chromosome_w_inter_disc_loops(:,2))*2 max(chromosome_w_inter_disc_loops(:,2))*2])
    zlim([min(chromosome_w_inter_disc_loops(:,3))*2 max(chromosome_w_inter_disc_loops(:,3))*2])
    set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
    xlabel('x','FontSize', 24)
    ylabel('y','FontSize', 24)
    zlabel('z','FontSize', 24)
    view(30,2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %output for CSynth: single cell contact probability
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PCSynth = zeros(nat_sum(size(P,1)),3);
    count = 0;
    for i=1:1:size(P,1)
        for j=i:1:size(P,1)
            count = count + 1;
            PCSynth(count,1)=i;
            PCSynth(count,2)=j;
            %PCSynth(count,1)=i*resolution;
            %PCSynth(count,2)=j*resolution;
            PCSynth(count,3)=P(i,j);
        end
    end

    PCSynth(:,3) = round(PCSynth(:,3),3,"significant");
    %%monomer coordinates in CSynth compatible format
    %writematrix(PCSynth,sprintf('cholesteric_CSynth_D4_%d.txt',d),'Delimiter','tab')
    %%monomer coordinates in xyz format
    %writematrix(chromosome_w_inter_and_intra_disc_loops,sprintf('cholesteric_monomer_locations_54_discs_%d.txt',d),'Delimiter','tab');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %figure S1B in paper
    %annotate location of cholesteric discs on contact map
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

    % f1=figure
    % %plot cholesteric disc track
    % hold on
    % histogram(discs_indicies, size(chromosome_w_inter_and_intra_disc_loops,1))
    % for i = 1:1:size(discs_ends,1)
    %     xline(discs_starts(i),'w','LineWidth',3)
    %     xline(discs_ends(i),'w','LineWidth',3)
    % end
    % title('Cholesteric Disc Locations')
    % box off
    % axis off
    % daspect([100 1 1]) %squish height of figure
    % hold off
    % ax1=gca;
    %
    % % f2=figure
    % % %plot intra-disc loops track
    % % histogram(intra_disc_loop_indicies, size(chromosome_w_inter_and_intra_disc_loops,1), 'EdgeColor', 'r')
    % % title('Intra Disc Loop Locations')
    % % box off
    % % axis off
    % % daspect([100 1 1]) %squish height of figure
    % % ax2=gca;
    %
    % f2=figure
    % %plot hic matrix
    % imagesc(P);
    % hold on;
    % colorbar
    % xlabel('primary sequence [bp]', 'fontsize', 24)
    % ylabel('primary sequence [bp]', 'fontsize', 24)
    % xlim([0 total_chromosome_length])
    % ylim([0 total_chromosome_length])
    % ax = gca;
    % ax.XTickLabel = ax.XTick*resolution;
    % ax.YTickLabel = ax.YTick*resolution;
    % set(gca,'ColorScale','log')
    % maj_axis_chr=colorbar;
    % maj_axis_chr.Label.String = 'Probability of Contact';
    % maj_axis_chr.FontSize = 18;
    % caxis([10^(-3) 10^0]);
    %
    % ax2=gca;
    %
    % figure
    % tcl=tiledlayout(8,1); % 8 by 1 tile layout
    % ax1.Parent=tcl;
    % ax1.Layout.Tile=1; % placing fig1 in the 1st tile
    % ax2.Parent=tcl;
    % ax2.Layout.Tile=2; % placing fig2 in the 2nd tile
    % % ax3.Parent=tcl;
    % % ax3.Layout.Tile=3; % placing fig3 in the 3nd tile
    % ax2.Layout.TileSpan=[7 1];  %fig2 spans tiles 2 through 8
    % linkaxes([ax1 ax2],'x')
    % xlim([0 size(chromosome_w_inter_and_intra_disc_loops,1)])
    % close(f1)
    % close(f2)

    %calculate contact probability curve
    %number of pairs of genomic positions separated by s on a given chromosome is Lc-s, where Lc is the length of the chromosome
    [s Ps] = contact_probability_xyz(chromosome_w_inter_and_intra_disc_loops,distance_threshold);

    s_avg{d} = s';
    Ps_avg{d} = Ps';

    %truncate end of chromosomes that are slightly longer, portions of
    %discs are not deleted
    if size(P,1)<total_chromosome_length
        P=padarray(P,[total_chromosome_length-size(P,1) total_chromosome_length-size(P,1)],0,'post');
    end

    P_agg = P_agg+P(1:total_chromosome_length,1:total_chromosome_length); %single cell contact maps are superimposed, then averaged
end

P_agg=P_agg./num_chroms;

%population level (aggregate) Hi-C contact map
figure
imagesc(P_agg);
hold on;
colorbar
xlabel('primary sequence [bp]', 'fontsize', 24)
ylabel('primary sequence [bp]', 'fontsize', 24)
xlim([0 total_chromosome_length])
ylim([0 total_chromosome_length])
ax = gca;
ax.XTickLabel = ax.XTick*resolution;
ax.YTickLabel = ax.YTick*resolution;
set(gca,'ColorScale','log')
maj_axis_chr=colorbar;
maj_axis_chr.Label.String = 'Contact Probability';
maj_axis_chr.FontSize = 18;
caxis([10^(-3) 10^0]);

s_agg=1:1:round(total_chromosome_length*0.95);
s_agg=s_agg';
Ps_agg = [];
for i=1:1:round(total_chromosome_length*0.95)
    Ps_agg = [Ps_agg; sum(diag(P_agg,i))/(total_chromosome_length-i)];
end

%shift y-axis so P(s=0) matches empirical dinoflagellate data, P(s=0)=10^-2
Ps_agg = Ps_agg./Ps_agg(1).*10^-2;

figure
hold on
ax=gca;
A=10;
xA=linspace(5*10^3,8*10^6,100);
plot(xA,A*xA.^(-0.5),'--k', 'Linewidth', 2)
plot(xA,A/100*xA.^(-0.2),'--r', 'Linewidth', 2)
set(ax,'xScale', 'log')
set(ax,'YScale', 'log')
plot(s_agg.*resolution,Ps_agg,'Color',[0.4 0.4 0.4],'Linewidth', 2)
xlabel('s','FontSize', 24)
ylabel('P(s)','FontSize', 24)
xline(2*(mean(mon_per_disc_agg)+mean(mean_loop_length_agg))*resolution,'--',{'2x Sequence',' in Layer'},'FontSize', 19)
xline(2*(mean(disc_diameter))*resolution,'--',{'2x Disc',' Diameter'},'FontSize', 19)
ax = gca;
ax.FontSize = 16;
ylim([10^-6 10^0])
xlim([10^4 2.5*10^7])

%output for CSynth: population level (aggregate) contact probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCSynth = zeros(nat_sum(round(total_chromosome_length*0.95)),3);
count = 0;
for i=1:1:round(total_chromosome_length*0.95)
    for j=i:1:round(total_chromosome_length*0.95)
        count = count + 1;
        %PCSynth(count,1)=i;
        %PCSynth(count,2)=j;
        PCSynth(count,1)=i*resolution;
        PCSynth(count,2)=j*resolution;
        PCSynth(count,3)=P_agg(i,j);
    end
end

PCSynth(:,3) = round(PCSynth(:,3),3,"significant");
%writematrix(PCSynth,strcat('cholesteric_short_frac_loops_0.1_CSynth_D4_aggregated_num_chroms_',num2str(num_chroms),'.txt'),'Delimiter','tab');

%calculate error bar for contact proability curve
%vertically concatenate single cell contact probability curves, pad with NaNs
maxNumCol = max(cellfun(@(c) size(c,2), Ps_avg));  % max number of columns
aMat = cell2mat(cellfun(@(c){padarray(c,[0,maxNumCol-size(c,2)],NaN,'Post')}, Ps_avg)');
colMeans = mean(aMat,1,'omitnan')';

Ps_avg = cellfun(@transpose,Ps_avg,'UniformOutput',false);
s_avg = cellfun(@transpose,s_avg,'UniformOutput',false);

[max_size, max_index] = max(cellfun('size', s_avg, 1));
[min_size, min_index] = min(cellfun('size', s_avg, 1));

% figure
% hold on
% ax=gca;
% err = [];
% for v = 1:1:size(s_avg{min_index},1)
%     spread_Ps = [];
%     for d = 1:1:num_chroms
%         spread_Ps = [spread_Ps Ps_avg{d}(v)];
%     end
%     err = [err std(spread_Ps)];
% end
% errorbar(s_avg{min_index}*5000,colMeans(1:size(s_avg{min_index},1)),err)
%legendStrings = "d_{cutoff} = " + string(distance_threshold) + "um";
%legend(legendStrings, 'Location', 'southwest')
% xlabel('Genomic separation [bp]', 'fontsize', 24)
% ylabel('Probability of contact', 'fontsize', 24)
% set(ax,'xScale', 'log')
% set(ax,'YScale', 'log')
% set(gca,'YLim',[10^-5 10^0],'YTick',10.^(-5:0))
% ax = gca;
% ax.FontSize = 22;
% xlim([10^4,10^7]);
% ylim([10^-5,10^0]);
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
