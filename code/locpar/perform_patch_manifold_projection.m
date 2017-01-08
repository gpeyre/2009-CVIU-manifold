function M1 = perform_patch_manifold_projection(M,D,options)

w = size(D,1);
m = size(D,3);
n = size(M,1);
 
D1 = reshape(D,[w*w m]);

H = compute_all_patch(M, w, options);
q = size(H,4); % number of patch
% remove mean
mu = mean(mean(H));
H = H - repmat(mu, [w w]);
% reshape
H1 = reshape(H,[w*w q]);

% compute inner product
A = H1' * D1;
% retrieve best correlation
[tmp,I] = max(abs(A),[],2);
C = A( (I-1)*q + (1:q)' );
% recompose
P = D(:,:,I) .* repmat(reshape( C, [1 1 q] ), [w w]);
P = reshape(P,[w w 1 q]);
% re-inject mean
P = P + repmat(mu, [w w]);
% reconstruct
M1 = compute_all_patch(P, n, options);