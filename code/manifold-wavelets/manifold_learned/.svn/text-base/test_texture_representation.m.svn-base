% test for the representation of a texture using another one

path(path, './nntools/');

rep = 'images/';
name = 'crochetcropti';
name = 'empilements';
name = 'frenchfries';

name0 = 'nuts';
name1 = 'mures';
name1 = 'frenchfries';


name0 = 'dunes';
name1 = 'nuts';

name1 = 'yellow_peppers';
name0 = 'tomatoes';

name0 = 'dunes';
name1 = 'dunes';

n = 64;
M0 = rescale( load_image([rep name0]) );
M1 = rescale( load_image([rep name1]) );
n0 = min(n, size(M0,1));
n1 = min(n, size(M1,1));
if strcmp(name0,name1)==0
    M0 = M0(end/2-n0/2+1:end/2+n0/2,end/2-n0/2+1:end/2+n0/2,:);
    M1 = M1(end/2-n1/2+1:end/2+n1/2,end/2-n1/2+1:end/2+n1/2,:);
    M00 = M0;
    % perform color equalization
    M0 = ycbcr2rgb( perform_histogram_matching(rgb2ycbcr(M0), rgb2ycbcr(M1), 255) );
    M0 = rescale(M0);
else
    %M0 = M0(1:n0,1:n0,:);
    M0 = M0(1:n0,end-n0+1:end,:);
    M1 = M1(end-n1+1:end,end-n1+1:end,:);
    M00 = M0;
end


w = 2;
options.w = w;
options.wmax = 12;
s = size(M1,3);

options.sampling = 'uniform';
% options.sampling = 'random';
% options.nbr = 4000;
options.use_nn = 1;
options.pca_numvecs = 50;
if 0
    [U,V,Err] = perform_texture_estimation(M0,M1,options);
else
    [U,V,Uevol,Vevol] = perform_texture_estimation_multiscale(M0,M1,options);
end

E = zeros(n1,n1,3);
E(:,:,1) = U;
E(:,:,2) = V;
E = rescale(E);

% perform image reconstruction
I = sub2ind(size(M0),U,V);
M1_rec = M1;
for i=1:size(M0,3)
    Mu = M0(:,:,i);
    M1_rec(:,:,i) = Mu(I);
end

clf;
subplot(2,3,1);
imagesc(M0); axis image; axis off;
title('Source texture');
subplot(2,3,2);
imagesc(M1); axis image; axis off;
title('Target texture');
subplot(2,3,3);
imagesc(M1_rec); axis image; axis off;
title('Reconstructed texture');
subplot(2,3,4);
imagesc(E); axis image; axis off;
title('Coordinate Mapping');
subplot(2,3,5);
if exist('Err')
    imagesc(Err); axis image; axis off;
    colormap gray(256);
    title('Patchwise error');
end

%%% save image
rep = ['results/texture_representation/' name0 '_' name1 '/'];
if ~exist(rep)
    mkdir(rep);
end
% save original image
warning off;
imwrite( M00, [rep 'exemplar_' name0 '_original.png'], 'png' );
imwrite( M0, [rep 'exemplar_' name0 '_equalized.png'], 'png' );
imwrite( M1, [rep 'target_' name1 '_original.png'], 'png' );
imwrite( M1_rec, [rep 'target_' name1 '_reconstructed.png'], 'png' );
imwrite( rescale(E), [rep 'mapping_' name1 '_to_' name0 '.png'], 'png' );
if exist('Err')
imwrite( rescale(Err), [rep 'error_' name1 '_to_' name0 '.png'], 'png' );
end
warning on;