clc
clear

%see SI: Wang, H., Yang, J., Zhang, Y., Qian, J., & Wang, J. (2022). Reconstruct high-resolution 3D genome structures for diverse cell-types using FLAMINGO. Nature Communications, 13(1), 2645.
coccoid = importdata('coccoid.xyz');
mastigote = importdata('mastigote.xyz');

bp = coccoid(:,1);
coccoid(:,1) = [];
mastigote(:,1) = [];

A=coccoid'*mastigote;
[U,S,V] = svd(A);
coccoid_aligned = coccoid*U*V';
writematrix([bp coccoid_aligned],'coccoid_aligned','Delimiter','tab')

% function [ rmse,A2 ] = SVD3D( A,B)
% 
% % Conformation alignment using the singular value decomposition (SVD) algorithm
% % References:
% % Sorkine,O. Least-squares rigid motion using SVD. (2009) Technical notes, 120, 3.
% n=size(A,1);
% [ret_R, ret_t] = RigidTransform3D(A, B);
% A2 = (ret_R*A') + repmat(ret_t, 1, n);
% A2 = A2';
% 
% % Find the error
% err = A2 - B;
% err = err .* err;
% err = sum(err(:));
% rmse = sqrt(err/n);
% 
% end