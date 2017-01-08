% test for some pca computation for patches

rep = 'images/';
name = 'dunes';
name = 'mures';
name = 'crochetcropti';
name = 'empilements';
name = 'frenchfries';
name = 'nuts';

M = load_image([rep name]);
n = min(128, size(M,1));
M = M(end/2-n/2+1:end/2+n/2,end/2-n/2+1:end/2+n/2,:);
w = 11;
s = size(M,3);
M = rescale(M);

rep = ['results/patches_pca/' name '/'];
if ~exist(rep)
    mkdir(rep);
end

% save original image
warning off;
imwrite( M, [rep name '_original.png'], 'png' );

% load random patches
options.sampling = 'uniform';
X = compute_patch_library(M,w,options);
p = size(X,4);
ww = 2*w+1;
Xs = reshape( X, [prod(size(X(:,:,:,1))) p] );

% pick the patch
if ~exist('x') || ~exist('y')
    clf;
    imagesc(M); axis image; axis off;
    [y,x,p] = ginput(1);
    x = round(x);
    y = round(y);
end
Mi = M(x-w:x+w,y-w:y+w,:);

% compute distance
D = compute_distance_to_points(Xs,Mi(:)); D = D(:);

% compute neihbbors
thresh = 0.03;
thresh = 0.05;
epsilon = max(D)*thresh;
nn = 80;
d = sort(D); 
epsilon = d(nn);
thresh = epsilon / max(D);

I = find(D<=epsilon);
DI = D(I);
Xsnn = Xs(:,I);
Xnn = X(:,:,:,I);
str = [name '_x' num2str(x) 'y' num2str(y) '_nn' num2str(round(100*thresh)) ];

disp(['Number of neighbors ' num2str(length(I)) ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display some typical patches in neighbors
clf;
typ = [0 0.2 0.5 0.7 1];
for i=1:length(typ)
    t = typ(i);
    [tmp,m] = min( abs(DI-epsilon*t) );
    A = Xnn(:,:,:,m);
    subplot(1,length(typ),i);
    imagesc(A); axis image; axis off;
    title([num2str(t) '%']);
    % save
    warning off;
    imwrite( A, [rep str '_neighb' num2str(typ(i)) '.png'], 'png' );
    warning on;
end


% compute pca vectors
numvecs = 20;
[Y,X1,v] = pca(Xsnn,numvecs);


warning off;
imwrite( Mi, [rep str '_patch' num2str(i) '.png'], 'png' );

nbr = 9;
figure; clf;
for i=1:nbr
    A = reshape( Y(:,i), [ww,ww,s] );
    subplot(3,3,i);
    imagesc(rescale(A));
    axis image; axis off;
    title(['Eigenvector ' num2str(i)]);
    % save image
    imwrite( rescale(A), [rep str '_pca' num2str(i) '.png'], 'png' );
end
warning on;

% display spectra
figure; clf;
nbr = 16;
plot( rescale( v(1:nbr) ), '.-' );
axis tight;
saveas(gcf, [rep str '_eigenvalues.eps'], 'epsc');