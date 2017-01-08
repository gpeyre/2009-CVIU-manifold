function [L,W] = compute_neighborhood_laplacian(M, sigma_sp, sigma_va)

if nargin<2
    sigma_sp = 5; % number of spacial pixels
end
if nargin<3
    M = rescale(M);
    sigma_va = 0.05; 
end

n = size(M,1);
[Y,X] = meshgrid(1:n,1:n);
W = zeros(n^2,n^2);
for i=1:n^2
    progressbar(i,n^2);
    d = (X-X(i)).^2 + (Y-Y(i)).^2;
    ws = exp( -d.^2/(2*sigma_sp^2) );
    d = (M-M(i)).^2 + (M-M(i)).^2;
    wv = exp( -d.^2/(2*sigma_va^2) );
    W(:,i) = ws(:) .* wv(:);
end
W = (W+W')/2;
d0 = sum(W,2);
d = sum(W,2); D = diag(d);
% symetrize
K = W ./ (sqrt(d*d'));
L = eye(n^2)-K;