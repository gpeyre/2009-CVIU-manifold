% test for 2D regular signal manifold

name = 'regular2d';

colmap = gray(256);
colmap = jet(256);

n = 50;
k = 3;
m = 25; % dimreduc
options.alpha = 3.5;

options.ndims = m;
options.k = k;
options.bound = 'per';

rep = ['results/manifold/' name '/'];
if not(exist(rep))
    mkdir(rep);
end
repeps = [rep '/eps/' ];
if not(exist(repeps))
    mkdir(repeps);
end

M = rescale( load_image('fnoise',n,options) );
nn = prod(size(M));

warning off;
imwrite(clamp(M), [rep name '-original.png'], 'png');
warning off;

H = perform_lowdim_embedding(M,options);
H = reshape(H, nn, m)';

%% compute 3D embedding
dx = M([2:end 1],:)-M([end 1:end-1],:);
dy = M(:,[2:end 1])-M(:,[end 1:end-1]);
dx = rescale(dx); dy = rescale(dy);
S = cat(3, dx,dy,M);
S(end+1,:,:) = S(1,:,:);
S(:,end+1,:) = S(:,1,:);

%% display embedding
clf;
col = M; col(end+1,:) = col(1,:); col(:,end+1) = col(:,1);
surf(S(:,:,1), S(:,:,2), S(:,:,3), col); 
view(3); axis equal;
camlight;
lighting phong;
xlabel('df(x)/dx');
ylabel('df(x)/dy');
zlabel('f(x)');
colormap(colmap);

saveas(gcf, [rep name '-manifold.png'], 'png');
saveas(gcf, [repeps name '-manifold.eps'], 'epsc');
saveas(gcf, [rep name '-nonlocal.png'], 'png');

%% display height field
clf;
x = linspace(0,1,n+1);
surf(x,x,S(:,:,3), col);
view(3); axis equal;
camlight;
lighting phong;
colormap(colmap);
saveas(gcf, [rep name '-semilocal.png'], 'png');

%% display 2D image
clf;
x = linspace(0,1,n+1);
surf(x,x,0*S(:,:,3), col);
axis([0 1 0 1 -0.1 0.1]); 
view(3); axis equal;
shading interp;
camlight;
lighting none;
colormap(colmap);
saveas(gcf, [rep name '-local.png'], 'png');

%% ask for some point
clf;
imagesc(M); axis image; axis off;
colormap(colmap);
title('Click on some salient feature');
[y,x] = ginput(1); x = round(x); y = round(y);
u = sub2ind(size(M),x,y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute some geodesic distances
options.nn_nbr = 4;
Gr = compute_nn_graph(H,options);
Gr(Gr==Inf) = 0; Gr = sparse(Gr);
% compute the distance between point_list and the remaining points
Dgeo = compute_distance_graph(Gr, u);
% compute the euclidean distance, which is the distance inside the big manifold
Deucl = sqrt( compute_distance_to_points(H,H(:,u)) );
Deucl = reshape(Deucl, n,n); Dgeo = reshape(Dgeo, n,n);

clf;
subplot(1,2,1);
imagesc(Deucl); axis image; axis off;
title('Eulidean distance');
subplot(1,2,2);
imagesc(Dgeo); axis image; axis off;
title('Geodesic distance');
colormap jet(256);
saveas(gcf, [rep name '-geodesic-dist-recap.png'], 'png');
warning off;
imwrite( apply_colormap(Dgeo,'jet'), [rep name '-geodesic-dist.png'], 'png' );
imwrite( apply_colormap(Deucl,'jet'), [rep name '-euclidean-dist.png'], 'png' );
warning on;