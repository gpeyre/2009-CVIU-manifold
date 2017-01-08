% test for building k-nn 

path(path, './nntools/');

rep = 'images/';
name = 'dunes';
name = 'mures';
name = 'crochetcropti';
name = 'empilements';
name = 'frenchfries';
name = 'nuts';

M = load_image([rep name]);
n = min(128, size(M,1));
M = M(end/2-n/2+1:end/2+n/2,end/2-n/2+1:end/2+n/2,:);
w = 12;
s = size(M,3);
M = rescale(M);

rep = ['results/patches_manifold/' name '/'];
if ~exist(rep)
    mkdir(rep);
end

% save original image
warning off;
imwrite( M, [rep name '_original.png'], 'png' );

% load random patches
options.sampling = 'uniform';
options.sampling = 'random';
options.nbr = 4000;
options.wmax = 4;
X = compute_patch_library(M,w,options);
p = size(X,4);
ww = 2*w+1;
Xs = reshape( X, [prod(size(X(:,:,:,1))) p] );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% precomputation to accelerate the nn search
% build k-nn
k = 11;
% perform pca simplification
numvecs = 30;
fprintf('--> Performing PCA reduction ... ');
[Y,Xs1,v,Psi] = pca(Xs,numvecs);
fprintf('done.\n');


fprintf('--> Pre-compute NN ... ');
% prepare for nn search
atria = nn_prepare(Xs1');
% perform search query
[index,distance] = nn_search(Xs1', atria, 1:p, k, 0);
fprintf('done.\n');