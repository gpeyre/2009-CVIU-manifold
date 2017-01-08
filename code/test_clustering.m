% test for clustering


name = 'brick';
name = 'square';
name = 'barb';
name = 'chessboard';
name = 'grass';
name = 'lines_orientation_smooth';
name = 'locally-parallel';
name = 'text';
name = 'corral';
name = 'reptilskin';
name = 'fingerprint';
name = 'fabric2';
name = 'dead_leaf';
name = 'locally-parallel';
name = 'locparx';
name = 'cross-dot-line';
name = 'filet';
name = 'cordes';
name = 'grid-circles';
name = 'patchwork';
path(path, '../images/');
path(path, 'images/');


n = 90;
nn = n^2;
k = 4;  % half size of the window

n0 = [];
mu = 0;
c = [];
switch name
    case 'fingerprint'
        n0 = 200;
    case 'chessboard'
        n0 = n;
        options.width = n/4;
        options.sigma = 1;
        mu = 0.003; % additional noise level
    case 'corral'
        n0 = 120;
        k=2;
    case 'dead_leaf'
        n0 = 350;
        k = 5;
    case 'fabric2'
        n0 = 200;
    case 'grass'
        n0 = 100;
    case 'lines_orientation_smooth'
        n0 = 80;
        k = 2;
    case 'reptilskin'
        n0 = 180;
        k = 5;
    case 'locally-parallel'
        k = 4;
    case 'text'
        k = 4;
        n0 = 300;
    case 'locparx'
        n0 = 200;
    case 'grid-circles'
        n0 = n;
        k = 5;
    case 'cross-dot-line'
        n0 = 100;
        k = 3;
    case 'filet'
        n0 = 220;
        c = [40 n0/2+5];
    case 'cordes'
        n0 = 120;
end
m = min((2*k+1)^2 - 4, 30); % for dimreduc

rep = ['results/clustering/' name '/'];
if not(exist(rep))
    mkdir(rep);
end
repeps = [rep '/eps/' ];
if not(exist(repeps))
    mkdir(repeps);
end

options.k = k;
kk = (2*k+1)^2;

options.width = 0.4;
options.frequency = 13;

%% load some shape
if not(strcmp(name, 'patchwork'))
    M = load_image(name,n0, options);
    M = sum(M,3);
    M = rescale( crop(M,n,c) );
else
    % names{end+1} = 'fine_grid';
    names{1} = 'reptilskin';
    names{2} = 'fine_grid';
    % names{3} = 'herringbone';
    names{3} = 'warped_grid';
    names{4} = 'hair';
    names{5} = 'fur';
    p = 4; Ma = {};
    for k=1:p
        n0 = [];
        switch names{k}
            case {'corral' 'dunes'}
                n0 = 200;
            case 'warped_grid'
                n0 = 350;
            case 'wood'
                n0 = 230;
            case 'grass-big'
                n0 = 200;
            case 'reptilskin'
                n0 = 180;
            case 'hair'
                n0 = 200;
        end
        A = load_image(names{k}, n0);
        A = sum(A,3);
        A = rescale( crop(A,n) );
        % perform histogram flattening
        A = perform_histogram_equalization( A, linspace(0,1,n^2) );
        Ma{k} = A;
    end
    [M,Id] = compute_texture_patchwork(Ma);
end

if strcmp(name, 'grid-circles')
    h = compute_gaussian_filter([11 11], 1.5/(2*n), [n n]);
    M = perform_convolution(M, h);
end

%% save original
clf;
imagesc(M); axis image; axis off;
colormap gray(256);
warning off;
imwrite(clamp(M), [rep name '-original.png'], 'png');
warning off;

%% compute patch vectors
options.ndims = m;
[H,Proj,Psi] = perform_lowdim_embedding(M,options);
H = reshape(H, nn, m)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute some geodesic distances
options.nn_nbr = 10;
Gr = compute_nn_graph(H,options);
Gr(Gr==Inf) = 0; Gr = sparse(Gr);

% compute the distance between point_list and the remaining points
Dgeo = compute_distance_graph(Gr, 1:nn);
% euclidean distance
Deucl = sqrt( compute_distance_matrix(H) );

options.niter = 8;
options.nrestart = 5;
nb = 4;
[Igeo, Vgeo] = perform_kmeans_distance(Dgeo,nb, options);
[Ieucl, Veucl] = perform_kmeans_distance(Deucl,nb, options);

Igeo = reshape(Igeo,[n n]);
Ieucl = reshape(Ieucl,[n n]);
clf;
subplot(1, 3, 1);
imagesc(M); axis image; axis off;
title('Image');
subplot(1, 3, 2);
imagesc(apply_colormap(Ieucl,'jet')); axis image; axis off;
title('Euclidean');
subplot(1, 3, 3);
imagesc(apply_colormap(Igeo,'jet')); axis image; axis off;
title('Geodesic');
colormap gray(256);


warning off;
imwrite( rescale(M), [rep name '-original.png'], 'png' );
imwrite( rescale(apply_colormap(Ieucl,'jet')), [rep name '-segm-eucl.png'], 'png' );
imwrite( rescale(apply_colormap(Igeo,'jet')), [rep name '-segm-geod.png'], 'png' );
warning on;

return;


