% test for geodesic strucuture of image patches
%
%   Copyright (c) 2007 Gabriel Peyre

% test for spectral graph computations
n = 40;

%% load image

path(path, '../images');

name = 'regular1d';
name = 'gaussiannoise';
name = 'step';
name = 'stepregular';

name = 'brick';
name = 'regular1';
name = 'square';
name = 'sparsecurves';
name = 'regular3';
name = 'locally-parallel';
name = 'diskregular';
name = 'barb';
name = 'fingerprint';
name = 'disk';

name = 'dead_leaf';
name = 'locally-parallel';
name = 'text';
name = 'lines_orientation_smooth';
name = 'grass';
name = 'chessboard';
name = 'fabric2';
name = 'potatoe';
name = 'grid-circles';
name = 'barb';
name = 'corral';
name = 'reptilskin';
name = 'disk';

path(path, 'toolbox/');
path(path, './images/');
path(path, '../images/');


options.sigma = 0;
options.alpha = 3;
n0 = n;
ndims = 2;
options.bound = 'sym';
name1 = name;
mu = 0; % additional noise level
switch name
    case {'disk' 'diskregular'}
        n0 = n;
        options.sigma = 1;
        mu = 0.03;
    case 'fingerprint'
        n0 = 200;
    case 'brick'
        n0 = 80;
    case 'chessboard'
        options.width = n/4;
        options.sigma = 1;
    case 'barb'
        n0 = [];
        mu = 0;
    case 'locally-parallel'
        n0 = [];
        mu = 0;
    case 'regular1d'
        ndims = 1;
        options.bound = 'per';
        name1 = 'regular';
        options.alpha = 1.5;
        mu = 0;
    case 'gaussiannoise'
        ndims = 1;
        options.bound = 'per';
        options.sigma = 256/10;
        mu = 0;
    case 'step'
        ndims = 1;
        options.bound = 'sym';
    case 'stepregular'
        mu = 0;
        ndims = 1;
        options.bound = 'sym';
    case 'corral'
        n0 = 120;
        n0 = 100;
    case 'dead_leaf'
        n0 = 350;
    case 'fabric2'
        n0 = 200;
    case 'grass'
        n0 = 100;
    case 'lines_orientation_smooth'
        n0 = 80;
        k = 2;
    case 'reptilskin'
        n0 = 180;
        % k = 5;
    case 'locally-parallel'
        k = 3;
    case 'text'
        k = 4;
        n0 = 300;
    case 'potatoe'
        k = 3;
        n0 = n;
end

rep = ['results/diffusion/' name '/' ];
if not(exist(rep))
    mkdir(rep);
end
repimg = [rep 'img/' ];
if not(exist(repimg))
    mkdir(repimg);
end

if ndims==2
    M = load_image(name,n0,options);
    M = sum(M,3);
    M = rescale( crop(M,n) );
    % add some noise to avoid degenerescence
else
    n = 256;
    M = load_signal(name1,n,options);
end
M = rescale(M+randn(size(M))*mu);

if strcmp(name, 'grid-circles')
    M = perform_convolution(M, ones(3)/9);
end


if ndims==1
    k = 4;
    m = 2*k+1;
else
    if not(exist('k'))
        k = 3;  % half size of the patches
    end
    m = 25; % dimreduc
end
nn = prod(size(M));


compute_noisy_laplacians = 1;
sigma = 0.06;
if compute_noisy_laplacians
    Anoise = randn(size(M))*sigma;
else
    Anoise = M*0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% displaying the image as a 3D surface
if ndims==2
    clf;
    subplot(1,2,1);
    imagesc(M); axis image; axis off;
    colormap gray(256);
    subplot(1,2,2);
    surf(M); axis tight;
    shading interp; lighting gouraud; camlight;
    saveas(gcf, [rep name '-image.png'], 'png');
    warning off;
    imwrite(clamp(M), [repimg name '-original.png'], 'png');
    warning off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute low dim embedding
disp('--> Computing low dimensional embedding');
options.ndims = m;
options.k = k;
if ndims==2
    H = perform_lowdim_embedding(M+Anoise,options);
    m = size(H,3);
    H = reshape(H, nn, m)';
else
    H = perform_lowdim_embedding_1d(M+Anoise,options);
end
% number of embeded points
q = min(nn,5000);
sel = randperm(nn); sel =  sel(1:q);
if q==nn
    sel = 1:q;
end
Hs = H(:,sel);
% central value mapping
Fval = M(sel);
% position mapping
if ndims==2
    [Y,X] = meshgrid(1:n,1:n);
else
    X = 1:n;
end
Fpos = X(sel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ask the user for some point
use_input = 0;
if ndims==2
    clf;
    imagesc(M); axis image; axis off;
    colormap gray(256);
    title('Click on some salient feature');
    if use_input
        [y,x] = ginput(1); x = round(x); y = round(y);
    else
        y = 10; x = 10;
    end
    u = sub2ind(size(M),x,y); u = sel(u);
    clf;
    imagesc(M); axis image; axis off;
    colormap gray(256);
    title('Click inside some object');
    if use_input
        [yin,xin] = ginput(1); xin = round(xin); yin = round(yin);
    else
        yin = 10; xin = 10;
    end
else
    clf;
    plot(M); axis tight;
    title('Click on some salient feature');
    [x,y] = ginput(1); x = round(x); y = round(y);
    u = x; u = sel(u);
    clf;
    plot(M); axis tight;
    title('Click inside some object');
    [xin,yin] = ginput(1); xin = round(xin); yin = round(yin);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute distance matrix
disp('--> Computing kernels');
sigma = 0.1;
Dist = compute_distance_matrix(Hs) / (2*k+1)^2;
Dist = (Dist+Dist')/2;   % enforce symmetry
% weights
W0 = exp( -Dist / sigma^2 );
d0 = sum(W0,2);
Llgd = {'un-normalized' 'normalized'};
Llgd = {'neigh' 'normalized' 'local'};
Llgd = {'neigh' 'un-normalized' 'normalized' 'local'};
Llgd = {'normalized'};
Llgd = {'neigh' 'un-normalized' 'local'};
Llgd = {'un-normalized'};
nlap = length(Llgd);
W = {}; L = {}; U = {}; S = {}; Us = {}; Ua = {}; D = {}; d = {}; K = {};
for i=1:nlap
    if strcmp(Llgd{i}, 'local')
        % local laplacian
        options.ndims = ndims;
        options.bound = 'per';
        L{i} = -compute_laplacian_matrix(n, options);
        L{i} = full(L{i});
        % display a row
        sigma_sp = 3;
        [Y,X] = meshgrid(1:n,1:n);
        dist = (X-x).^2 + (Y-y).^2;
        w = exp( -dist./(sigma_sp^2) );
        clf;
        imagesc(w); axis image; axis off;
        colormap gray(256);
        warning off;
        imwrite( rescale(w), [repimg name '-' Llgd{i} '-weights.png'], 'png');
        warning on;
    else
        if strcmp(Llgd{i}, 'neigh')
            % neighboorhood
            sigma_sp = 4.5;
            sigma_va = 0.1;
            [tmp,W{i}] = compute_neighborhood_laplacian(M+Anoise, sigma_sp, sigma_va);
        else
            W{i} = W0;
        end
        if strcmp(Llgd{i}, 'normalized')
            % normalize to density
            W{i} = W{i} ./ (d0*d0');
        end
        d{i} = sum(W{i},2); D{i} = diag(d{i});
        % display a row
        w = reshape( W{i}(u,:), n,n);
        w(u) = max( w(w~=w(u)) );
        clf;
        imagesc(rescale(w).^0.33); axis image; axis off;
        colormap gray(256);
        warning off;
        imwrite( rescale(w), [repimg name '-' Llgd{i} '-weights.png'], 'png');
        warning on;
        % symetrize
        K{i} = W{i} ./ (sqrt(d{i}*d{i}'));
        L{i} = eye(q)-K{i};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spectral decomposition
disp('--> Performing spectral decomposition');
for i=1:nlap
    disp(['* ' Llgd{i} '.']);
    [U{i},S{i}] = eig(L{i});
    S{i} = diag(S{i});
    Us{i} = U{i}; Ua{i} = U{i}';
    if not(strcmp(Llgd{i}, 'local'))
        % synthesis mapping u^->u = D^{-1/2}*U
        Us{i} = Us{i} ./ repmat( sqrt(d{i}(:)), [1 q] );
        % analysis mapping u->u^ = U'*D^{1/2}
        Ua{i} = Ua{i} .* repmat( sqrt(d{i}(:)'), [q 1] );
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% display some eigenvectors
a = 6; b = 8;
for k=1:nlap
    if ndims==2
    if strcmp(Llgd{k}, 'local')
        ilist = 8 + (0:a*b-1);
    else
        ilist = round( linspace(2,350,a*b) );
    end
    else
        ilist = 2:20;
    end
    clf;
    for i=1:a*b
        v = U{k}(:,ilist(i));
        if ndims==2
            v = reshape( v, n,n);
            if strcmp(Llgd{k}, 'local')
                w = rescale(v);
            else
                w = clamp(0.4*v/std(v(:)),-1,1);
                w = (w+1)/2;
            end
        else
            w = v;
        end
        subplot(a,b,i);
        if ndims==2
            image( w*255 ); axis image; axis off;
            warning off;
            imwrite( rescale(v), [repimg name '-' Llgd{k} '-eigv-' num2str(i) '.png'], 'png');
            warning on;
        else
            plot(w); axis tight;
        end
        title(['Eigv.' num2str(ilist(i))]);
    end
    colormap gray(256);
    saveas(gcf, [rep name '-eigenvectors-' Llgd{k} '.png'], 'png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute some diffusion
tlist = {};
tlist{1} = [1 2 4 8 16 32]*0.5;   % diffusion times
tlist{2} = [0 1 2 4 8 16]*0.3;    % denoising times
switch name
    case 'diskregular'
%        tlist = [1 2 4 8 16 32]*16;
    case 'barb'
        tlist{1} = [1 2 4 8 16 32]*4;
%        tlist{2} = [0 1 2 4 8 16]*0.5;
    case 'chessboard'
        % tlist{1} = [1 2 4 8 16 32]*64; 
    case 'fingerprint'
        tlist{1} = [1 2 4 8 16 32]*12;
    case 'locally-parallel'
        tlist{1} = [1 2 4 8 16 32]*12;
end

% initial conditions
if ndims==2
    Mdirac = zeros(n); Mdirac(xin,yin) = 1;
else
    Mdirac = zeros(n,1); Mdirac(xin) = 1;
end
% noisy image
eta = 0.05;
sigma = 0.08; mu = 0.05;
if ndims==2
    Mn = rescale(M,mu,1-mu) + randn(n)*sigma;
else
    Mn = rescale(M,mu,1-mu) + randn(n,1)*sigma;
end
% intial conditions
f = {Mdirac, Mn};
test = {'diffusion', 'denoising'};

for a=1:length(f)
for k=1:nlap
    disp(['--> Computing diffusion ' Llgd{k}]);
    % compute eigendecomposition
    fa = f{a}(:);
    fu = Ua{k}*fa;
    % plot the coefficients
    clf;
    plot( clamp( abs(U{k}'*fa), 0,1.5) ); axis tight;
    if strcmp(test{a}, 'denoising')
        saveas(gcf, [repimg name '-' test{a} '-' Llgd{k} '-coefficients.png'], 'png');
    end
    clf;
    for i=1:min(length(tlist{1}),6)
        % compute diffusion
        t = tlist{a}(i);
        w = exp(-t*S{k}); % diffusion weights
        fi = Us{k} * ( fu .* w );
        subplot(2,3,i);
        if ndims==2
            fi = reshape(fi,n,n);
            if a==1
                A = clamp(fi); % apply_colormap(clamp(fi),'jet');
            else
                A = rescale(fi);
            end
            imagesc(fi); axis image; axis off;
            colormap gray(256);
            warning off;
            imwrite(rescale(A), [repimg name '-' test{a} '-' Llgd{k} '-' num2str(i) '.png'], 'png');
            warning off;
        else
            plot(fi); axis tight;
        end
        title(['t=' num2str(t)]);
    end
    saveas(gcf, [rep name '-' test{a} '-' Llgd{k} '.png'], 'png');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare sparse approximation in spectral and wavelets
Mw = [];
lgd = {};
for i=1:nlap
    Mw(:,end+1) = U{i}'*M(:);
    lgd{end+1} = Llgd{i};
end

% enlage the image a little bit to mach 
if 0
nn = 2^ceil(log2(n));
MM = [M; M(end:-1:2*end-nn+1,:)];
MM = [MM MM(:,end:-1:2*end-nn+1)];
MM = M;

Jmin = 0;
options.wavelet_vm = 4;
options.wavelet_type = 'daubechies';
MW = perform_wavelet_transform(MM, Jmin, +1, options);
MW = MW(randperm(nn^2)); MW = MW(1:n^2);
MW = MW * sqrt( sum(M(:).^2) / sum(MW(:).^2) );
Mw(:,end+1) = MW(:);
lgd{end+1} = 'wavelets';
end

MW2 = cumsum( sort( Mw.^2 ) );
MW2 = MW2(end:-1:1,:);
sel = ceil(0.01*n^2):round(0.7*n^2);
clf;
t = {'k:', 'b', 'r--', 'g-.'};
hold on;
for i=1:size(MW2,2)
    h = plot(100*sel/n^2, log2(MW2(sel,i)), t{i}); axis tight;
    set(h, 'LineWidth',3);
end
hold off;
% legend(lgd);
saveas(gcf, [rep name '-approx-err-decay.png'], 'png');
saveas(gcf, [rep name '-approx-err-decay.eps'], 'epsc');

