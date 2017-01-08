function [U,V,Uevol,Vevol] = perform_texture_estimation_multiscale(M0,M1,options)

% perform_texture_estimation - fit the patches of a given texture onto another one.
%
%   [U,V,Err] = perform_texture_estimation(M0,M1,options);
%
%   M0 and M1 are images 
%
%   You can specify the size of the patches in options.w.
%   If you want to use PCA dimension reduction for speeding up the
%       computation, then set options.pca_numvecs to the number of reduced
%       dimension (default 15). Set to -1 to use full dimension (no speed up).
%   For using fast closest point queries, you have to installe nntools/ 
%       in your path. To override the use of nntools, you can set
%       options.use_nn=0
%
%   Copyright (c) Gabriel Peyré 2006


options.null = 0;
if isfield(options, 'w0')
    w0 = options.w0;
else
    w0 = 1;
end
if isfield(options, 'Jmax')
    Jmax = options.Jmax;
else
    Jmax = 4;
end
if isfield(options, 'k0')
    k0 = options.k0;
else
    k0 = 4;
end
if isfield(options, 'pca_numvecs')
    m = options.pca_numvecs;
else
    m = 15;
end

% perform initial estimation
options.w = 2^Jmax*w0;
options.wmax = 2^2*w0;
options.sampling = 'uniform';
[U,V,Err] = perform_texture_estimation(M0,M1,options);
Uevol{1} = U; Vevol{1} = V;

n1 = size(M1,1);

for j=Jmax-1:-1:0
    w = w0*2^j;
    k = k0*2^j;
    % load patches library
    [X0,X,Y] = compute_patch_library(M0,w,options);
    X1 = compute_patch_library(M1,w,options);
    p0 = size(X0,4); p1 = size(X1,4);
    X0s = reshape( X0, [prod(size(X0(:,:,:,1))), p0] );
    X1s = reshape( X1, [prod(size(X1(:,:,:,1))), p1] );
    % perform dimension reduction
    m = min(m,size(X0s,1));
    [P,tmp,v,Psi] = pca(X0s,m);
    X0s = X0s - repmat(Psi,[1,size(X0s,2)]);
    X1s = X1s - repmat(Psi,[1,size(X1s,2)]);
    X0s = P' * X0s;
    X1s = P' * X1s;
    % estimate the k-nn of current points 
    Ind = sub2ind(size(M0),U,V);
    % prepare for nn search
    fprintf('--> Pre-compute NN ... ');
    atria = nn_prepare(X0s');
    [I,distance] = nn_search(X0s', atria, Ind(:), k, 0);
    fprintf('done.\n');
    % compute the distances
    for i=1:n1^2
        % find closest point to X1s(:,i) in set X0s(:,I(i,:))
        Ind = I(i,:);
        D = compute_distance_to_points( X0s(:,Ind), X1s(:,i) );
        [tmp,J] = min(D);
        U(i) = X( Ind(J) );
        V(i) = Y( Ind(J) );
    end
    Uevol{end+1} = U; Vevol{end+1} = V;
end