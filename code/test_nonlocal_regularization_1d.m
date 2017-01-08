% test for non-local regularization
n = 1024;
name = 'piece-regular';

rep = 'results/nonlocal-1d/';
if not(exist(rep))
    mkdir(rep);
end
repeps = [rep 'eps/'];
if not(exist(repeps))
    mkdir(repeps);
end

f = rescale( load_signal(name, n) );

% extract patches
w1 = 3; w = 2*w1+1;
dx = -w1:w1;
I = repmat(dx,[n 1]) + repmat((1:n)',[1 w]);
I(I<=0) = 2-I(I<=0);
I(I>n) = 2*n-I(I>n);

H = f(I');

% compute euclidian distance
D = sqrt( compute_distance_matrix(H)/w );

sigma = 0.05;
W = exp( -(D/sigma).^2 );

warning off;
imwrite(rescale(W), [name  '-weights.png'], 'png');
imwrite(rescale(D), [name '-dist.png'], 'png');
warning on;


d = diag(sum(W,2));
% renormalized smoothing operator
W1 = W ./ repmat( sum(W,2), [1 n]);
% laplacian operator
L = d-W; L = (L+L')/2;
% normalized laplacian
L1 = eye(n)-W1;


%% compute eigendecomposition
disp('--> Computing eigendecomposition');
tic;
[V,delta] = eig(L); delta=real(diag(delta));
toc;


%% compute hard and soft thresholdings
fW = V'*f;
T = 0.15;
fWT = perform_thresholding(V'*f, T, 'hard');

plot(fW); axis([1 n -.5 .5]);
saveas(gcf, [rep name '-transformed.png'], 'png');
saveas(gcf, [repeps name '-transformed.eps'], 'epsc');
plot(fWT); axis([1 n -.5 .5]);
saveas(gcf, [rep name '-thresholded.png'], 'png');
saveas(gcf, [repeps name '-thresholded.eps'], 'epsc');


%% test for thresholding flow
Tlist = linspace(0,1,3); % [0 0.2 0.2 0.3 0.4];
Fhard = []; Fsoft = []; Fquad = []; Fgauss = [];
for i=1:length(Tlist)
    t = Tlist(i);
    Fhard(:,end+1) = V*perform_thresholding(V'*f, t, 'hard');
    Fsoft(:,end+1) = V*perform_thresholding(V'*f, t, 'soft');
    A = diag( 1./(1+delta*t) );
    Fquad(:,end+1) = V*A*V'*f;
    A = diag( exp( -delta*t ) );
    Fgauss(:,end+1) = V*A*V'*f;
end
my_plot(Fhard, 2, 'hard');
my_plot(Fsoft, 2, 'soft');
my_plot(Fquad, 2, 'quad');
my_plot(Fgauss, 2, 'gauss');