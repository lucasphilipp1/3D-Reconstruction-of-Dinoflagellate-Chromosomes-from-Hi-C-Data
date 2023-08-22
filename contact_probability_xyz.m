function [s Ps] = contact_probability_xyz(model,distance_threshold)
Ps=zeros(size(model,1)-1,1);
s=1:1:size(model,1)-1;
s=s';
D = pdist(model);
D = squareform(D);
cont_count_mat = D;
cont_count_mat(cont_count_mat<distance_threshold) = 0;
cont_count_mat(cont_count_mat>distance_threshold) = 1;
cont_count_mat(:) = ~cont_count_mat;
cont_count=zeros(1,size(model,1)-1);
for i=1:1:size(model,1)-1
    cont_count(i) = sum(diag(cont_count_mat,i));
    Ps(i)=cont_count(i)/(size(model,1)-i);
end
end