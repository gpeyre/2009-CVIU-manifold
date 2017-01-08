% test for 1D step edges
name = 'step1d';

rep = ['results/manifold/' name '/' ];
if not(exist(rep))
    mkdir(rep);
end
repeps = [rep '/eps/' ];
if not(exist(repeps))
    mkdir(repeps);
end

n = 1024;

%% step function with given heights
h = [2 3 5 1];
h = repmat(h,[n/length(h) 1]);
f = rescale(h(:));

%% set of parameters (manifold)
k = 17;
nc = 20;
na = 20;
nb = 20;
clist = linspace(-k+1,k-1,nc);
alist = linspace(0,1,na);
blist = linspace(0,1,nb);
[A,B,C] = ndgrid(alist,blist,clist);
I = find(A+B/2<=1 & A-B/2>=0);
A = A(I); B = B(I); C = C(I);
m = length(A); % number of basis signals

%% set of patches
x = (-k:k)';
H = repmat(x, [1 m]);
H = double(H>=repmat(C',[2*k+1 1]))-1/2;
H = H .* repmat(B', [2*k+1 1]) + repmat(A', [2*k+1 1]);
% remove doublons
H = unique(H', 'rows')';

%% dimension reduction on H
ndims = 3;
options.nn_nbr = 7;
xy = isomap(H,ndims,options);

%% compute the path
options.bound = 'sym';
options.k = k;
Hf = perform_lowdim_embedding_1d(f,options);
d = compute_distance_to_points(H,Hf);
[tmp,nn_list] = min(d');
nn_list(diff(nn_list)==0)=[];


%% display
clf;
hold on;
plot_scattered(xy, H(k+1,:));
h = plot3(xy(nn_list,1), xy(nn_list,2), xy(nn_list,3), 'k');
hold off;
set(h, 'LineWidth', 4);
view(3);
saveas(gcf, [rep name '-manifold.png'], 'png');
saveas(gcf, [repeps name '-manifold.eps'], 'epsc');