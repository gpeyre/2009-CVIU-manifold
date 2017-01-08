function [U,V,Err] = perform_texture_estimation(M0,M1,options)

% perform_texture_estimation - fit the patches of a given texture onto another one.
%
%   [U,V,Err] = perform_texture_estimation(M0,M1,options);
% or 
%   [U,V,Err] = perform_texture_estimation(X0,M1,options);
% or 
%   UV = perform_texture_estimation(X0,M1,options);
%
%   M0 and M1 are images 
%   X0 is a set of p patches of size (2*w+1,2*w+1,s) where s is the number of colors, 
%       so X0 should be of size (2*w+1,2*w+1,s,p).
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

if isfield(options, 'use_nn')
    use_nn = options.use_nn;
else
    if exist('nn_prepare')
        use_nn = 1;
    else
        use_nn = 0;
        warning('nntools is missing, you should add them in the path.');
    end
end

options.null = 0;
if isfield(options, 'w')
    w = options.w;
else
    w = 7;
end

if size(M0,4)>1
    % use provided patches
    X0 = M0;
else
    % extract from image
    [X0,X,Y] = compute_patch_library(M0,w,options);
end

n1 = size(M1,1);
options.sampling = 'uniform';
X1 = compute_patch_library(M1,w,options);

p0 = size(X0,4); p1 = size(X1,4);
X0s = reshape( X0, [prod(size(X0(:,:,:,1))), p0] );
X1s = reshape( X1, [prod(size(X1(:,:,:,1))), p1] );
d = size(X0s,1);

% perform dimension reduction
if isfield(options, 'pca_numvecs')
    m = options.pca_numvecs;
else
    m = 15;
end
m = min(m,d);
if m>0 && m<d
    fprintf('--> Performing dimension reduction ...');
    q = min(5000,size(X0s,2));
    sel = randperm(size(X0s,2));  sel = sel(1:q);
    [P,tmp,v,Psi] = pca(X0s(:,sel),m);
    X0s = X0s - repmat(Psi,[1,size(X0s,2)]);
    X1s = X1s - repmat(Psi,[1,size(X1s,2)]);
    X0s = P' * X0s;
    X1s = P' * X1s;
    fprintf('done.\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use fast NN Search
if use_nn
    k = 11;
    fprintf('--> Pre-compute NN ... ');
    % prepare for nn search
    atria = nn_prepare(X0s');
    fprintf('done.\n');
    fprintf('--> Performing estimation ...');
    % perform search query
    [I,distance] = nn_search(X0s', atria, X1s', 1, 0);
    U = reshape(X(I),n1,n1);
    V = reshape(Y(I),n1,n1);
    Err = reshape(distance,n1,n1);
    fprintf('done.\n');
    if nargout==1
        UV = zeros([size(U),2]);
        UV(:,:,1) = U; UV(:,:,2) = V;
        U = UV;
    end
    return;
end

U = zeros(n1);
V = zeros(n1);
Err = zeros(n1);
fprintf('--> Performing estimation ...');
h = waitbar(0,'Computing texture estimation ...');
for i=1:n1^2
    waitbar(i/n1^2, h);
    % compute distance to patch
    if use_nn
        % use fast query
    else
        % use slow matlab
        D = compute_distance_to_points(X0s,X1s(:,i));
        [tmp,I] = min(D);
    end
    U(i) = X(I);
    V(i) = Y(I);
    Err(i) = sqrt( tmp );
end
close(h);
fprintf('done.\n');

if nargout==1
    UV = zeros([size(U),2]);
    UV(:,:,1) = U;
    UV(:,:,2) = V;
end