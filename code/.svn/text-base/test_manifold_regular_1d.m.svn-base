% test for 1D regular signal manifold

name = 'regular1d';

rep = ['results/manifold/' name '/' ];
if not(exist(rep))
    mkdir(rep);
end
repeps = [rep '/eps/' ];
if not(exist(repeps))
    mkdir(repeps);
end

n = 1024;
options.bound = 'per';
options.sigma = 100 * n/512;


f = load_signal('gaussiannoise',n,options);
f = rescale(f(:));
x = linspace(0,1,n);

% compute derivatives
df = f([2:end 1])-f([end 1:end-1]);
df = rescale(df);
h = [f, df]';

% plot the function
clf;
hh = plot(x, f, 'k');
set(hh, 'LineWidth', 3);
axis tight;
saveas(gcf, [rep name '-function.png'], 'png');
saveas(gcf, [repeps name '-function.eps'], 'epsc');
% plot the curve
clf;
hold on;
scatter(f,df,50,x,'filled');
hh = plot(f,df,'k');
set(hh, 'LineWidth', 2);
plot([0 1 1 0 0], [0 0 1 1 0], 'k');
hold off;
axis([0 1 0 1]);
axis square;
axis tight;
colormap jet(256);
saveas(gcf, [rep name '-manifold.png'], 'png');
saveas(gcf, [repeps name '-manifold.eps'], 'epsc');

%% ask for some point
clf;
plot(f); axis tight;
title('Click to select some point');
[xin,yin] = ginput(1); xin = round(xin); yin = round(yin);

%% compute the geodesic distance
k = 3;
options.bound = 'per';
options.k = k;
H = perform_lowdim_embedding_1d(f,options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute some geodesic distances
options.nn_nbr = 4;
Gr = compute_nn_graph(H,options);
Gr(Gr==Inf) = 0; Gr = sparse(Gr);
% compute the distance between point_list and the remaining points
Dgeo = compute_distance_graph(Gr, xin);
% compute the euclidean distance, which is the distance inside the big manifold
Deucl = sqrt( compute_distance_to_points(H,H(:,xin)) );

clf;
hh = plot(x,Dgeo, 'b', x, Deucl, 'r--'); axis tight;
set(hh, 'LineWidth', 3);
legend('Euclidean', 'Geodesic');
saveas(gcf, [rep name '-geodesic-dist.png'], 'png');
saveas(gcf, [repeps name '-geodesic-dist.eps'], 'epsc');


