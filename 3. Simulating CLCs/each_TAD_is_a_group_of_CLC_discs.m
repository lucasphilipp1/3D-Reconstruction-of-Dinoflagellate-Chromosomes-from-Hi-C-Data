%clc
%clear
%Author: Lucas Philipp
%iteratively create CLC chromosomes
%discs order is permuted randomly within groups, recapitulating TADs
%simulate CLC Hi-C contact maps/contact probability curves
%
%requires: statistics and machine learning toolbox

%s_agg_original = s_agg;
%Ps_agg_original = Ps_agg;

num_chroms = 350-190-48-106; %number of simulated chromosomes
resolution = 5000; %number of base pairs per monomer

%contact probability curve data
s_avg = cell(1, num_chroms); %cell array for genomic separation
Ps_avg = cell(1, num_chroms); %cell array for contact probability

total_chromosome_length = 4000; %number of monomers
%P_agg = zeros(total_chromosome_length); %collect Hi-C contact maps for different chromosomes

disc_diameter = []; %collect disc_diameter for different chromosomes

TAD_bp = [
    21449      3136855; %TAD1 start/stop
    3138645    5456633; %TAD2 start/stop
    5440909    7815782; %etc...
    7817867    8660723;
    8677664    10342262;
    10355695   14413730;
    14436643   19268194
    ]; %TAD start/stop positions in base pairs

total_chromosome_length_bp = 19282064;

TAD_mon_index = round(TAD_bp./total_chromosome_length_bp*total_chromosome_length); %convert TAD positions from bp to monomer index

xy_spacing = 0.0627; %arbitrary spatial distance unit, same as iter_Cholesteric_HiC

for d = 1:1:num_chroms
    d
    chromosome=[]; %position of DNA
    %7.655Mbp/(10 layers)∗1monomer/5kbp≈(150 monomers)/layer
    num_discs = 30; %target number of cholesteric discs

    %at least 15% of total monomers shold be inter-disc loops in order for
    %connections to be long enough (as connected discs are no longer
    %immeadiately above or below)
    frac_tot_sequence_in_loops = 0.2; %Fraction
    frac_loop_sequence_inter_disc = 0.75; %Fraction

    % Find the target total length of sequence in core and in loops
    total_inter_disc_loop_length = round((total_chromosome_length)*frac_tot_sequence_in_loops*frac_loop_sequence_inter_disc);
    total_intra_disc_loop_length = round((total_chromosome_length)*frac_tot_sequence_in_loops*(1-frac_loop_sequence_inter_disc));

    num_mon_per_disc = round((total_chromosome_length - total_inter_disc_loop_length - total_intra_disc_loop_length)/num_discs);

    num_mon_per_disc_stack = round((TAD_mon_index(:,2)-TAD_mon_index(:,1))*(1-frac_tot_sequence_in_loops)); %all CLC monomers not in loops are in discs

    %cholesteric pitch P: length in microns along long axis of the chromosome
    %corresponding to a full turn of the nematic director
    pitch = 1.30/2; %in microns

    min_axis_chr=0.5; %radius of cholesteric disc, in microns

    % optional slanted discs (discs not perpendicular to chromosome long axis)
    %tilt_angle = 10; %in degrees, 0 is perpendicular to chromosome long axis
    %tilt_angle = random('Normal',0,30); %in degrees
    tilt_angle = 0;

    %add noise to z coordinate from normal dist with 0 mean and std dev z_spacing
    noise_z = false;

    xy_spacing = min_axis_chr*sqrt(pi/num_mon_per_disc); %mesh spacing in microns
    z_spacing = xy_spacing;

    %distance cutoff to count as HiC contact
    distance_threshold = xy_spacing*2;

    %2D grid of points
    x = -1:xy_spacing:1; %mesh in microns
    y = -1:xy_spacing:1; %mesh in microns

    dtheta_layer=360*z_spacing/pitch; %DNA fibres in a layer are rotated dtheta_layer degrees counterclockwise relative to the previous layer

    disc_groupings = round(num_mon_per_disc_stack/num_mon_per_disc);
    chol_layers = sum(disc_groupings);
    maj_axis_chr=chol_layers*z_spacing; %in microns

    z = linspace(-chol_layers/2*z_spacing,chol_layers/2*z_spacing,chol_layers);

    layer_order = 1:chol_layers;

    shuffled_layer_order = [];
    index = 1;

    % Process each group
    for i = 1:length(disc_groupings)
        % Extract group
        group = layer_order(index : index + disc_groupings(i) - 1);

        % Randomize within group
        group = group(randperm(disc_groupings(i)));

        % Append to output
        shuffled_layer_order = [shuffled_layer_order, group];

        % Update index
        index = index + disc_groupings(i);
    end

    %discs are different sizes, need cell array
    discs = {};

    intra_disc_loop_potential_starts_ind = [];
    for i = 1:1:chol_layers

        [X,Y] = meshgrid(x,y);

        r = min_axis_chr; %discs have constant width throughout chromosome (no taper at chromosome ends)

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
        R=[cosd((shuffled_layer_order(i)-1)*dtheta_layer) -sind((shuffled_layer_order(i)-1)*dtheta_layer); sind((shuffled_layer_order(i)-1)*dtheta_layer) cosd((shuffled_layer_order(i)-1)*dtheta_layer)]; %create 2D rotation matix, applies to (X,Y) pairs
        rotcoord=[];
        for j=1:1:size(xEllipse,1)
            rotcoord=[rotcoord; [xEllipse(j), yEllipse(j)]*R']; %multiply by rotation matrix
        end

        %create array of potential start positions for intra-disc extra chromosomal loops
        intra_disc_loop_potential_starts_ind = [intra_disc_loop_potential_starts_ind; size(chromosome,1) + flip_these];

        chromosome = [chromosome; [rotcoord(:,1), rotcoord(:,2), z(shuffled_layer_order(i)).*ones(size(xEllipse))]]; %append positions of DNA fibres for new cholesteric disc
        discs{end+1} = [rotcoord(:,1), rotcoord(:,2), z(shuffled_layer_order(i)).*ones(size(xEllipse))];
    end

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
    % zlim([min(chromosome(:,3))*2 max(chromosome(:,3))*2])
    % set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
    % for i = 0:0.05:1
    %     zlim([max(chromosome(:,3))*(i-0.03) max(chromosome(:,3))*(i+0.03)])
    %     pause(5)
    % end
    % view(0,90)
    % hold off

    %INTER DISC LOOPS
    idx=find(abs(diff(chromosome(:,3)))>0); %new disc starts when z coordinate changes. now next disc can be above or below previous disc
    num_inter_disc_loops = size(idx,1);
    mean_loop_length = total_inter_disc_loop_length/num_inter_disc_loops;
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

    %re-order inter_disc_loop_lengths so longer inter disc loops can connect discs located further away in z
    %%%
    order_discs = abs(diff(shuffled_layer_order))+rand(size(abs(diff(shuffled_layer_order)))).*1e-10; %breaks ties, but doesn't change order
    order_discs = (tiedrank(-order_discs));

    order_loop_lengths = inter_disc_loop_lengths+rand(size(inter_disc_loop_lengths)).*1e-10; %breaks ties, but doesn't change order
    order_loop_lengths = (tiedrank(-order_loop_lengths));

    [~, sort_idx] = ismember(order_discs, order_loop_lengths);
    inter_disc_loop_lengths = inter_disc_loop_lengths(sort_idx);
    %%%

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
        %reject walks with any steps of length spacing*5 or greater
        j=1;
        while j<=size(inter_disc_loops{i},1)-1
            while sqrt((inter_disc_loops{i}{j,1}-inter_disc_loops{i}{j+1,1}).^2 + (inter_disc_loops{i}{j,2}-inter_disc_loops{i}{j+1,2}).^2 + (inter_disc_loops{i}{j,3}-inter_disc_loops{i}{j+1,3}).^2) > z_spacing*5
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
    while sum_inter_disc_loops >= total_inter_disc_loop_length
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

    while sum(intra_disc_loop_lengths) <= total_intra_disc_loop_length
        new_loop = round(exprnd(mean_loop_length));
        % Reject loops lengths of 0
        while new_loop == 0
            new_loop = round(exprnd(mean_loop_length));
        end
        if sum(intra_disc_loop_lengths)+new_loop > total_intra_disc_loop_length
            break;
        end
        intra_disc_loop_lengths = [intra_disc_loop_lengths new_loop];
    end

    if sum(intra_disc_loop_lengths) < total_intra_disc_loop_length
        intra_disc_loop_lengths = [intra_disc_loop_lengths total_intra_disc_loop_length-sum(intra_disc_loop_lengths)];
    end

    intra_disc_loop_lengths = round(intra_disc_loop_lengths);

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
        %reject walks with any steps of length spacing*5 or greater
        j=1;
        while j<=size(intra_disc_loops{i},1)-1
            while sqrt((intra_disc_loops{i}{j,1}-intra_disc_loops{i}{j+1,1}).^2 + (intra_disc_loops{i}{j,2}-intra_disc_loops{i}{j+1,2}).^2 + (intra_disc_loops{i}{j,3}-intra_disc_loops{i}{j+1,3}).^2) > z_spacing*5
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

    chromosome_w_inter_and_intra_disc_loops_before_tilt_and_noise = chromosome_w_inter_and_intra_disc_loops;
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
    if noise_z == true
        chromosome_w_inter_and_intra_disc_loops(:,3) = chromosome_w_inter_and_intra_disc_loops(:,3) + normrnd(0,z_spacing,size(chromosome_w_inter_and_intra_disc_loops,1), 1);
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
    % view(90,10)
    % shuffled_layer_order'

    %
    % for i=skip+1:skip:numPoints
    %     clf(f)
    %     plot3(chromosome_w_inter_and_intra_disc_loops(:,1),chromosome_w_inter_and_intra_disc_loops(:,2),chromosome_w_inter_and_intra_disc_loops(:,3),'Color', [.6 .6 .6])
    %     xlim([min(chromosome_w_inter_and_intra_disc_loops(:,1))*1.1 max(chromosome_w_inter_and_intra_disc_loops(:,1))*1.1])
    %     ylim([min(chromosome_w_inter_and_intra_disc_loops(:,2))*1.1 max(chromosome_w_inter_and_intra_disc_loops(:,2))*1.1])
    %     zlim([min(chromosome_w_inter_and_intra_disc_loops(:,3))*1.1 max(chromosome_w_inter_and_intra_disc_loops(:,3))*1.1])
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
    %writematrix(chromosome_w_inter_and_intra_disc_loops,sprintf('CLC_each_TAD_is_a_group_of_discs_%d.txt',d),'Delimiter','tab');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %figure S1B in paper
    %annotate location of cholesteric discs on contact map
    all_disc_locations = [];
    intra_disc_loops_locations = [];
    inter_disc_loops_locations = [];

    for i = 1:numel(discs)
        all_disc_locations = [all_disc_locations; discs{i}];
    end

    idx = ismember(chromosome_w_inter_and_intra_disc_loops_before_tilt_and_noise, all_disc_locations, 'rows');
    disc_indices = find(idx);

    for i = 1:numel(inter_disc_loops)
        inter_disc_loops_locations = [inter_disc_loops_locations; cell2mat(inter_disc_loops{i})];
    end

    idx = ismember(chromosome_w_inter_and_intra_disc_loops_before_tilt_and_noise, inter_disc_loops_locations, 'rows');
    inter_disc_loops_indices = find(idx);

    for i = 1:numel(intra_disc_loops)
        intra_disc_loops_locations = [intra_disc_loops_locations; cell2mat(intra_disc_loops{i})];
    end

    idx = ismember(chromosome_w_inter_and_intra_disc_loops_before_tilt_and_noise, intra_disc_loops_locations, 'rows');
    intra_disc_loops_indices = find(idx);

    %calculate contact probability curve
    %number of pairs of genomic positions separated by s on a given chromosome is Lc-s, where Lc is the length of the chromosome
    [s Ps] = contact_probability_xyz(chromosome_w_inter_and_intra_disc_loops,distance_threshold);

    s_avg{d} = s';
    Ps_avg{d} = Ps';

    %old way of aligning chromosomes, start of first disc is always aligned
    %%%
    %truncate end of chromosomes that are slightly longer, portions of
    %discs are not deleted
    % if size(P,1)<total_chromosome_length
    %     P=padarray(P,[total_chromosome_length-size(P,1) total_chromosome_length-size(P,1)],0,'post');
    % end

    %P_agg = P_agg+P(1:total_chromosome_length,1:total_chromosome_length); %single cell contact maps are superimposed, then averaged
    %%%

    %new way of aligning chromosomes, if matrix is larger/shorter a random
    %crop/pad is applied to the left side, the crop/pad on the right side
    %is determined by the target chromosome size
    %%%
    size_diff = size(P,1) - total_chromosome_length;
    r = randi(abs(size_diff)+1)-1;
    if size_diff<0
        P=padarray(P,[r r],0,'pre');
        P=padarray(P,[abs(size_diff)-r abs(size_diff)-r],0,'post');
    end

    if size_diff>0
        P=P(r+1:size(P,1)-(size_diff-r),r+1:size(P,1)-(size_diff-r));
    end

    P_agg = P_agg + P;
end

P_agg=P_agg./num_chroms;

bar_data_indices = [];
label_size = 10;
TAD_boundary_index = [label_size+1; TAD_mon_index(1:end-1,2); size(P_agg,1)-(label_size+1)];
for i = 1:size(TAD_boundary_index, 1)
    startVal = TAD_boundary_index(i, 1)-label_size;
    endVal = TAD_boundary_index(i, 1)+label_size;

    % Generate the linspace row and convert it to a column
    linspaceVec = (startVal:1:endVal)';

    % Append to result
    bar_data_indices = [bar_data_indices; linspaceVec];
end

bar_data = zeros(1, size(P_agg,1));
bar_data(bar_data_indices) = -1;
bar_thickness = 200; %in pixels
bar_matrix = repmat(bar_data, bar_thickness, 1);
annotation_and_CM=[bar_matrix; P_agg];

figure;
imagesc(annotation_and_CM);
set(gca, 'ColorScale', 'log');
colormap(parula);
caxis([1e-3 1e0]);

hold on;

%can use two different colormaps for the same imagesc(); plot
% mask CLC disc track
mask_white = annotation_and_CM == 0;
mask_black = annotation_and_CM == -1;

% Overlay white
[rows, cols] = find(mask_white);
scatter(cols, rows, 1, 'w', 'filled');

% Overlay black
[rows, cols] = find(mask_black);
scatter(cols, rows, 1, 'k', 'filled');
xlabel('primary sequence [bp]', 'fontsize', 24)
ylabel('primary sequence [bp]', 'fontsize', 24)
xlim([0 total_chromosome_length])
ylim([0 total_chromosome_length])
ax = gca;
set(gca,'ColorScale','log')
maj_axis_chr=colorbar;
maj_axis_chr.Label.String = 'Contact Probability';
maj_axis_chr.FontSize = 18;
caxis([10^(-3) 10^0]);

ax.XTickLabel = ax.XTick*resolution;
%offset because of CLC disc track
ax.YTick = ax.YTick + bar_thickness;
ax.YTickLabel = (ax.YTick-bar_thickness)*resolution;
ytick_labels = string(ax.YTickLabel);

%starting ytick should be 0
idx = find(yticks == bar_thickness);
if ~isempty(idx)
    ytick_labels(idx) = "0";
end
ax.YTickLabel = ytick_labels;
ylim([0 size(annotation_and_CM,1)])

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
plot(xA,A/3*xA.^(-0.5),'--k', 'Linewidth', 2)
plot(xA,A/235*xA.^(-0.2),'--r', 'Linewidth', 2)
set(ax,'xScale', 'log')
set(ax,'YScale', 'log')
plot(s_agg.*resolution,Ps_agg,'Color',[0.4 0.4 0.4],'Linewidth', 2)
xlabel('s','FontSize', 24)
ylabel('P(s)','FontSize', 24)
xline(2*(mean(num_mon_per_disc)+mean(mean_loop_length))*resolution,'--',{'2x Sequence',' in Layer'},'FontSize', 19)
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
        PCSynth(count,1)=i;
        PCSynth(count,2)=j;
        %PCSynth(count,1)=i*resolution;
        %PCSynth(count,2)=j*resolution;
        PCSynth(count,3)=P_agg(i,j);
    end
end

PCSynth(:,3) = round(PCSynth(:,3),3,"significant");
%writematrix(PCSynth,strcat('cholesteric_short_frac_loops_0.1_CSynth_D4_aggregated_num_chroms_',num2str(num_chroms),'.txt'),'Delimiter','tab');
writematrix(PCSynth,strcat('CLC_permuted_fanc.txt'),'Delimiter','tab');

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
% ylabel('Contact Probability', 'fontsize', 24)
% set(ax,'xScale', 'log')
% set(ax,'YScale', 'log')
% set(gca,'YLim',[10^-5 10^0],'YTick',10.^(-5:0))
% ax = gca;
% ax.FontSize = 22;
% xlim([10^4,10^7]);
% ylim([10^-5,10^0]);
% grid on

writematrix(P_agg, 'CLC_permuted_contact_map.txt', 'Delimiter', 'tab')

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
