% test for texture modification


% size of input
n0 = 128;
n1 = 128*2;

rep = '';
rep = 'images/';
names = {};
% names{end+1} = 'fine_grid';
names{1} = 'reptil_skin';
names{2} = 'fine_grid';
% names{3} = 'herringbone';
names{3} = 'warped_grid';
names{4} = 'hair';
names{5} = 'fur';

name0 = names{1};
name0 = 'grass';
name0 = 'flow2';
name0 = 'reptil_skin';
names = {};
names{1} = 'simoncelli2';
names{1} = 'corral';
names{1} = 'simoncelli11';

repimg = 'results/modifications/';
if ~exist(repimg)
    mkdir(repimg);
end
repeps = 'results/modifications/eps/';
if ~exist(repeps)
    mkdir(repeps);
end

p = length(names);
M0 = {}; M1 = {};
for k=1:p
    name = names{k};
    if strcmp(name, 'corral') || strcmp(name, 'dunes')
        A = load_image([rep name], 256);
    elseif strcmp(name, 'warped_grid');
        A = load_image([rep name], 450);
    elseif strcmp(name, 'wood');
        A = load_image([rep name], 230);
    elseif strcmp(name, 'grass-big');
        A = load_image([rep name], 200);
    else
        A = load_image([rep name]);
    end
    A = sum(A,3);
    A = rescale(A);
    M0{k} = A(1:n0,1:n0);
    M1{k} = A(end-n1+1:end,end-n1+1:end);
    % M1{k} = A(1:n1,1:n1);
end

[M1a,Id] = compute_texture_patchwork(M1);

% M0a = M0{1};
A = load_image([rep name0]);
M0a = A(1:n0,1:n0);
M0a = sum( rescale( M0a ), 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform estimation
options.sampling = 'uniform';
options.use_nn = 1;
options.pca_numvecs = 12;
options.w = 2;
UV = perform_texture_estimation(M0a,M1a,options);

E = zeros(n1,n1,3);
E(:,:,1:2) = rescale(UV);

M1_rec = perform_texture_reconstruction(M0a,UV);

clf;
subplot(2,2,1);
imagesc(M0a); axis image; axis off;
title('Source texture');
subplot(2,2,2);
imagesc(M1a); axis image; axis off;
title('Target texture');
subplot(2,2,3);
imagesc(M1_rec); axis image; axis off;
title('Reconstructed texture');
subplot(2,2,4);
imagesc(E); axis image; axis off;
title('Coordinate Mapping');
colormap gray(256);

%%% perform texture analysis %%%
options.dim = 4;
Jmin = 2;
D = perform_texture_transform(M0a,UV,Jmin,+1,options);

if 0
UVb = perform_texture_transform(M0a,D,Jmin,-1,options);
M1_rec0 = perform_texture_reconstruction(M0a,UVb);
E = zeros(n1,n1,3);
E(:,:,1:2) = rescale(UVb);
end

clf;
plot_texture_wavelet(D, Jmin, options);
saveas(gcf, [repeps 'wavelets_0.eps'], 'epsc');

Ds = sqrt(sum(D.^2,3));

% perform various thresholding
T_list = [0 0.05 0.1 0.15 0.2 0.3 0.7 1 2 3 4 10];
warning off;
for i=1:length(T_list)
    % threshold
    t = T_list(i);
    I = find(Ds<t);
    D1 = D; 
    for s=1:size(D1,3)
        a = D1(:,:,s); a(I) = 0;
        D1(:,:,s) = a;
    end
    D1(1:2^Jmin,1:2^Jmin,:) = D(1:2^Jmin,1:2^Jmin,:);
    % bwd transform
    UVb = perform_texture_transform(M0a,D1,Jmin,-1,options);
    M1_rec = perform_texture_reconstruction(M0a,UVb);
    E = zeros(n1,n1,3); E(:,:,1:2) = rescale(UVb);
    
    % display
    if 0
    subplot(2,2,1);
    imagesc(M1_rec0); axis image; axis off;
    title('Original Image');
    subplot(2,2,2);
    imagesc(M1_rec); axis image; axis off;
    title('Modified');
    subplot(2,2,3);
    plot_texture_wavelet(D1, Jmin, options);
    subplot(2,2,4);
    imagesc(E); axis image; axis off;
    title('Coordinate Mapping');
    colormap gray(256);
    end
    
    % save
    imwrite(rescale(M1_rec), [repimg 'reconstruction' num2str(i) '.png'], 'png');
    imwrite(rescale(E), [repimg 'mapping' num2str(i) '.png'], 'png');
    clf;
    plot_texture_wavelet(D1, Jmin, options);
    saveas(gcf, [repeps 'wavelets_' num2str(i) '.eps'], 'epsc');
end

imwrite(rescale(M0a), [repimg 'im-src.png'], 'png');
imwrite(rescale(M1a), [repimg 'im-tgt.png'], 'png');

warning on;
