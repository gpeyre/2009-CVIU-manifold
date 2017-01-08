% test for extraction of eigenvectors of laplacian

n = 40;
n0 = n;

name = 'disk';

options.sigma = 3; % little blur


k = 4;  % half size of patches
q = 2*k+1;
nn = prod(size(M));

    
M = load_image(name,n0,options);
M = rescale( crop(M,n) );
M = M + randn(n)*.001;  % jitter a little

% compute low dimensional features
options.ndims = min(q^2,30);  % dimension reduction
options.k = k;
H = perform_lowdim_embedding(M,options);
m = size(H,3);
H = reshape(H, nn, m)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute distance matrix
sigma = 0.1;
Dist = compute_distance_matrix(H) / (2*k+1)^2;
Dist = (Dist+Dist')/2;   % enforce symmetry
% weights
W = exp( -Dist / (2*sigma^2) );
d0 = sum(W,2);
L = diag(d0) - W;

tic;
[U,S] = eig(L);
toc;
S = diag(S);

%tic;
%[U,S,V] = svd(L);
%toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display some eigenvectors
a = 6; b = 8;
ilist = round( linspace(2,350,a*b) );
clf;
for i=1:a*b
    v = U(:,ilist(i));
    v = reshape( v, n,n);
    w = clamp(0.4*v/std(v(:)),-1,1);
%    w = (w+1)/2;
    imageplot( v, ['Eigv.' num2str(ilist(i))], a,b,i );
    subplot(a,b,i);
end
