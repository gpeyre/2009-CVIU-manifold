% test for 2D step edges
name = 'step2d';

n = 80;
m = 30; % for dimreduc
options.n_delta = 30;
options.n_theta = 48;
options.sigma = 0.1;
nn = n^2;

rep = ['results/manifold/' name '/'];
if not(exist(rep))
    mkdir(rep);
end
repeps = [rep '/eps/' ];
if not(exist(repeps))
    mkdir(repeps);
end

k = 5;
options.k = k;
kk = (2*k+1)^2;
options.rescale = 0;
[P,delta_list,theta_list] = compute_edge_patches(2*k+1, options);
P = (P+1)/2;
HP = reshape(P, kk, size(P,3));

%% load some shape
M = load_image('images/potatoe', n);
warning off;
imwrite(clamp(M), [rep name '-original.png'], 'png');
warning off;

%% compute the mapping
options.ndims = m;
[H,Proj,Psi] = perform_lowdim_embedding(M,options);
HM = reshape(H, nn, m)';

%% compute best match 
HP = HP - repmat( Psi, [1 size(HP,2)] );
HP = Proj'*HP;
D = compute_distance_to_points(HP,HM);
[tmp,nn_list] = min(D, [], 2);
nn_list = reshape(nn_list,n,n);
Theta = theta_list( nn_list );
Delta = delta_list( nn_list );

Theta(Delta==min(Delta(:))) = 0;
Theta(Delta==max(Delta(:))) = pi;

%% display 
clf;
subplot(1,2,1);
imagesc(Delta); axis image; axis off;
title('Delta');
subplot(1,2,2);
imagesc(Theta); axis image; axis off;
title('Theta');
colormap jet(256);
saveas(gcf, [rep name '-theta-delta.png'], 'png');
warning off;
imwrite( apply_colormap(Theta,'jet'), [rep name '-theta.png'], 'png' );
imwrite( apply_colormap(Delta,'jet'), [rep name '-delta.png'], 'png' );
warning on;

%% dimension reduction on H
ndims = 3;
options.nn_nbr = 7;
xy = isomap(HP,ndims,options);

%% compute 1D path
nnpath = nn_list(60,:);
nnpath(diff(nnpath)==0)=[];

%% display
col = theta_list;
col = P(k+1,k+1,:); col = col(:);
clf;
hold on;
plot_scattered(xy, col);
h = plot3(xy(nnpath,1), xy(nnpath,2), xy(nnpath,3), 'k');
hold off;
set(h, 'LineWidth', 4);
view(3); axis tight;
saveas(gcf, [rep name '-manifold.png'], 'png');
saveas(gcf, [repeps name '-manifold.eps'], 'epsc');


%% ask the user for some point
n = 60;
M = load_image('images/potatoe', n);
[H,Proj,Psi] = perform_lowdim_embedding(M,options);
HM = reshape(H, n^2, m)';

clf;
imagesc(M); axis image; axis off;
colormap gray(256);
title('Click on some feature');
[y,x] = ginput(1); x = round(x); y = round(y);
u = sub2ind(size(M),x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute some geodesic distances
options.nn_nbr = 10;
Gr = compute_nn_graph(HM,options);
Gr(Gr==Inf) = 0; Gr = sparse(Gr);
% compute the distance between point_list and the remaining points
Dgeo = compute_distance_graph(Gr, u);
Dgeo = reshape(Dgeo,n,n);

clf;
imagesc(Dgeo); axis image; axis off;
colormap jet(256);
title('Geodesic distance');
colormap jet(256);
saveas(gcf, [rep name '-geodesic-dist-recap.png'], 'png');
warning off;
imwrite( apply_colormap(Dgeo,'jet'), [rep name '-geodesic-dist.png'], 'png' );
warning on;

