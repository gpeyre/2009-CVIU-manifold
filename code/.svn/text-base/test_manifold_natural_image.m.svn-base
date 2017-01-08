% test for natural image

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
name = 'grid-circles';
path(path, '../images/');
path(path, 'images/');


n = 60;
nn = n^2;
k = 4;  % half size of the window

n0 = [];
mu = 0;
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
        
end
m = min((2*k+1)^2 - 4, 30); % for dimreduc

rep = ['results/manifold/' name '/'];
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
options.frequency = 15;

%% load some shape
M = load_image(name,n0, options);
M = rescale( crop(M,n) );
M = rescale(M+randn(size(M))*mu);

if strcmp(name, 'grid-circles')
    h = compute_gaussian_filter([11 11], 1.5/(2*n), [n n]);
    M = perform_convolution(M, h);
end

%% save original
clf;
imageplot(M);
warning off;
imwrite(clamp(M), [rep name '-original.png'], 'png');
warning off;

%% compute the mapping
options.ndims = m;
[H,Proj,Psi] = perform_lowdim_embedding(M,options);
H = reshape(H, nn, m)';

%% dimension reduction on H
options.nn_nbr = 7;
xy = isomap(H,3,options);


%% 3D display
col = M(:);
clf;
plot_scattered(xy, col);
colormap jet(256);
view(3); axis tight;
saveas(gcf, [rep name '-manifold-3d.png'], 'png');
saveas(gcf, [repeps name '-manifold-3d.eps'], 'epsc');

%% 2D display
options.sampling = 'uniform';
P = compute_patch_library(M,k,options);
clf;
plot_flattened_dataset(xy(:,1:2),squeeze(P),25);
colormap gray(256);
saveas(gcf, [rep name '-manifold-2d.png'], 'png');
saveas(gcf, [repeps name '-manifold-2d.eps'], 'epsc');


%% ask the user for some point
clf;
imagesc(M); axis image; axis off;
colormap gray(256);
title('Click on some feature');
[y,x] = ginput(1); x = round(x); y = round(y);
u = sub2ind(size(M),x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute some geodesic distances
options.nn_nbr = 10;
Gr = compute_nn_graph(H,options);
Gr(Gr==Inf) = 0; Gr = sparse(Gr);
% compute the distance between point_list and the remaining points
Dgeo = compute_distance_graph(Gr, u);
Dgeo = reshape(Dgeo,n,n);
% euclidean distance
Deucl = sqrt( compute_distance_to_points(H,H(:,u)) );
Deucl = reshape(Deucl,n,n);

clf;
subplot(1,2,1);
imagesc(Dgeo); axis image; axis off;
title('Geodesic distance');
subplot(1,2,2);
imagesc(Deucl); axis image; axis off;
title('Euclidean distance');
colormap jet(256);
saveas(gcf, [rep name '-geodesic-dist-recap.png'], 'png');

warning off;
imwrite( apply_colormap(Dgeo,'jet'), [rep name '-geodesic-dist.png'], 'png' );
imwrite( apply_colormap(Deucl,'jet'), [rep name '-euclidean-dist.png'], 'png' );
warning on;

