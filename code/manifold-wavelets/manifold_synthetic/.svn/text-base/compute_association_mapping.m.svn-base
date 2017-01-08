function [A,D] = compute_association_mapping(H)

% compute_association_mapping - compute an association mapping
%
%   [A,D] = compute_association_mapping(H);
%
%   H is a (w,w,p) or a (w,w,s,p) matrix where 
%       w is the size of the patch, s is the number of color 
%       p is the number of patch.
%
%   k=A(i,j) describes that H(:,:,i) and H(:,:,j), when
%       mixed together in an horizontal direction, 
%       H(:,:,k).
%   D(i,j) states how well this mapping succeded.
%
%   Copyright (c) 2006 Gabriel Peyré

w = size(H,1);
w1 = (w-1)/2;
if size(H,4)>1
    % color mode
    s = size(H,3);
    p = size(H,4);
else
    % Gray level mode
    p = size(H,3);
    s = 1;
    H = reshape( H, [w w 1 p] );
end


% perform dimension reduction 
m = 15;
Hs = reshape(H,s*w^2,p);
[P,tmp,v,Psi] = pca(Hs,m);
Hs = Hs - repmat(Psi,[1,size(Hs,2)]);
Hs = P' * Hs;

% perform mixing
pmax = 100^2; % to avoid memory overflow
iter = ceil( p^2/pmax );
[J,I] = meshgrid(1:p,1:p); J = J(:); I = I(:);
A = []; D = [];
for i=1:iter
    
    sel = (i-1)*pmax+1:min(i*pmax,p^2);
    pp = length(sel);
    
    HH = zeros( w,w, s, pp );
    HH(1:w1+1,:,:,:) = HH(1:w1+1,:,:,:) + H(w1+1:end,:,:,I(sel));
    HH(w1+1:end,:,:,:) = HH(w1+1:end,:,:,:) + H(1:w1+1,:,:,J(sel));
    % fix middle column
    HH(w1+1,:,:,:) = HH(w1+1,:,:,:)/2;

    % perform dimension reduction
    HHs = reshape(HH,s*w^2,pp);
    HHs = HHs - repmat(Psi,[1,size(HHs,2)]);
    HHs = P' * HHs;

    % perform matching
    options.Y = HHs;
    [D1,A1] = compute_nn_distance(Hs, 1, options);
    A = [A;A1(:)];
    D = [D;D1(:)];

end

A = reshape(A,p,p);
D = reshape(D,p,p);