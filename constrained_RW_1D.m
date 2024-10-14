function walk = constrained_RW_1D(start,finish,requested_steps,step_size)
%constrainted random walk in one dimension
%walk does not include start point

%actual steps might be +1 requested steps
%e.g. 3 steps, end-start = 2, step_size 1, is impossible
%R and L need to be integers

R = (requested_steps + (finish-start)/step_size)/2;
L = (requested_steps - (finish-start)/step_size)/2;

Rint = abs(R-round(R)) < step_size/1e4;
Lint = abs(L-round(L)) < step_size/1e4;

if Rint==false | Lint==false
    steps=requested_steps+1;
else
    steps=requested_steps;
end

%R is number of right moves +1*step_size
%L is number of left moves -1*step_size
%R+L = steps
%step_size*(R-L) = end-start
%solving this system:

R = (steps + (finish-start)/step_size)/2;
L = (steps - (finish-start)/step_size)/2;

neg = -1*ones(round(L),1);
pos = +1*ones(round(R),1);
moves = [neg; pos];
moves=moves(randperm(size(moves,1)));

walk=NaN*zeros(steps,1);
current_pos=start;
for i=1:steps
    current_pos = current_pos + step_size*moves(i);
    walk(i)=current_pos;
end
end
