function walk = constrained_self_avoiding_RW_3D(start,finish,steps,step_size)
%contrained self avoiding 3D random walk
isunique=0;
while isunique==0
    x = constrained_RW_1D(start(:,1),finish(:,1),steps,step_size);
    y = constrained_RW_1D(start(:,2),finish(:,2),steps,step_size);
    z = constrained_RW_1D(start(:,3),finish(:,3),steps,step_size);

    size_x=size(x,1);
    size_y=size(y,1);
    size_z=size(z,1);

    %some paths require and even number of steps along a certain axis
    %some paths require an odd number of steps along a certain axis
    if find(find([size_x; size_y; size_z] == max([size_x; size_y; size_z]))==1) & size_x > steps
        x(end)=[];
    end
    if find(find([size_x; size_y; size_z] == max([size_x; size_y; size_z]))==2) & size_y > steps
        y(end)=[];
    end
    if find(find([size_x; size_y; size_z] == max([size_x; size_y; size_z]))==3) & size_z > steps
        z(end)=[];
    end
    walk = [x y z];

    %nearly self avoiding 3D, some overlap allowed, otherwise simulation time is long
    if steps > 100
        if size(unique(walk, "rows"),1) >= steps*0.75
            isunique=1;
        end
    elseif steps > 50
        if size(unique(walk, "rows"),1) >= steps*0.80
            isunique=1;
        end
    elseif steps > 10
        if size(unique(walk, "rows"),1) >= steps*0.85
            isunique=1;
        end
    else
        if size(unique(walk, "rows"),1) >= steps-1
            isunique=1;
        end
    end
end
end