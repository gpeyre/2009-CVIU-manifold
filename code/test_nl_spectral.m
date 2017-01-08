% test for spectral graph computations

n = 40;
m = 18; % dimreduc
k = 3;  % half size of the patches
nn = n^2;

%% load image
n0 = [];

repimg = '';
repimg = '../images/';

path(path, '../images');

name = 'barb';
name = 'fingerprint';
name = 'brick';
name = 'disk';

switch name
    case 'disk'
        n0 = n;
    case 'fingerprint'
        n0 = 200;
    case 'brick')
        n0 = 80;
end

test = 'cosine';
test = 'nlmeans';

repbase = ['results/eigendecomposition/' name '/' ];
rep = [repbase test '/'];
if not(exist(rep))
    mkdir(rep);
end

M = load_image(name,n0);
M = sum(M,3);
if strcmp(name, 'disk')
    h = compute_gaussian_filter( [9 9], 1/(2*n), [n n]);
    M = perform_convolution(M,h);
end
M = rescale( crop(M,n) );

if strcmp(test, 'nlmeans')

    %% compute low dim embedding
    options.ndims = m;
    options.k = k;
    H = perform_lowdim_embedding(M,options);
    H = reshape(H, nn, m);
    disp('--> Computing distance matrix');
    D = compute_distance_matrix(H') / (2*k+1)^2;
    % compute weights
    sigma = 0.05;
    W = exp( -D/sigma^2 );
    d = diag(sum(W,2));
    % renormalized smoothing operator
    W1 = W ./ repmat( sum(W,2), [1 nn]);
    % laplacian operator
    L = d-W; L = (L+L')/2;
    % normalized laplacian
    L1 = eye(nn)-W1;

else
    
    % classical laplacian
    L = -compute_laplacian_matrix(n, 'per');
    L = full(L);
    L1 = L;

end

%% compute eigendecomposition
disp('--> Computing eigendecomposition');
tic;
[V,delta] = eig(L); delta=real(diag(delta));
toc;

%% display some eigenvectors
warning off;
imwrite( rescale(M), [rep name '-original.png'], 'png');
warning on;

a = 3; b = 4;
if strcmp(test, 'nlmeans')
    ilist = 40 + (0:a*b-1);
else
    ilist = 4 + (0:a*b-1);
end

clf;
for i=1:a*b
    subplot(a,b,i);
    v = reshape( V(:,ilist(i)), n,n); 
    if strcmp(test, 'nlmeans')
        w = clamp(0.4*v/std(v(:)),-1,1);
        w = (w+1)/2;
    else
        w = rescale(v);
    end
    image( w*255 ); axis image; axis off;
    warning off;
    imwrite( rescale(v), [rep name '-eigv-' num2str(i) '.png'], 'png');
    warning on;
end
colormap gray(256);
saveas(gcf, [repbase name '-eigenvectors-' test '.png'], 'png');

% make some noise
sigma = 0.05;
mu = 0.05;
f = rescale(M,mu,1-mu) + randn(n)*sigma;

%% perform smooting
slist = linspace(1,0.05,12);
img = {}; coefs = {};
for i=1:length(slist)
    % project the function to approximate
    pf = V'*f(:);
    % set to zeros the high frequencies
    q = round(slist(i)*nn); % number of frequencies to keep
    pf(q+1:end) = 0;
    % reconstruct
    img{end+1} = clamp( reshape( V*pf, n,n ) );
    coefs{end+1} = clamp( abs(pf), 0,1.5 );
    clf;
    plot( coefs{i} ); axis tight;
    saveas(gcf, [rep name '-denoising-coefs-' num2str(i) '.png'], 'png');
end

% display images
clf;
for i=1:length(slist)
    subplot(3,4,i);
    image( img{i}*255); axis image; axis off;
    warning off;
    imwrite( img{i}, [rep name '-denoised-' num2str(i) '.png'], 'png');
    warning on;
end
colormap gray(256);
saveas(gcf, [repbase name '-denoising-imgs-' test '.png'], 'png');
% display coefs
clf;
for i=1:length(slist)
    subplot(3,4,i);
    plot( coefs{i} ); axis tight;
    warning off;
    imwrite( img{i}, [rep name '-denoised-' num2str(i) '.png'], 'png');
    warning on;
end
colormap gray(256);
saveas(gcf, [repbase name '-denoising-coefs-' test '.png'], 'png');

