% test for manifold approximation using local pca

path(path,'nntools/');

rep = 'images/';
name = 'mures';
name = 'empilements';
name = 'frenchfries';
name = 'nuts';
name = 'crochetcropti';
name = 'dunes';

M = load_image([rep name]);
n = min(128, size(M,1));
M = M(end/2-n/2+1:end/2+n/2,end/2-n/2+1:end/2+n/2,:);
w = 11;
s = size(M,3);
M = rescale(M);


w = 2;
% load random patches
options.sampling = 'uniform';
options.nbr = 2000;
X = compute_patch_library(M,w,options);
p = size(X,4);
ww = 2*w+1;
Xs = reshape( X, [prod(size(X(:,:,:,1))) p] );
% perform dimension reduction
m = 25;
[P,Xs,v,Psi] = pca(Xs,m);

[LocP,LocPsi,Segm] = perform_local_pca(Xs, options);