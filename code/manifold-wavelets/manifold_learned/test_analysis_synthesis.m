% test for texture analysis and synthesis
% test for the representation of a texture using another one

path(path, './nntools/');
rep = 'images/';

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


w = 4;
options.w = w;
options.wmax = 4;
s = size(M1,3);

options.sampling = 'uniform';
options.use_nn = 1;
options.pca_numvecs = 50;
UV = perform_texture_estimation(M0,M1,options);

E = zeros(n1,n1,3);
E(:,:,1:2) = rescale(UV);

M1_rec = perform_texture_reconstruction(M0,UV);

clf;
subplot(2,2,1);
imagesc(M0); axis image; axis off;
title('Source texture');
subplot(2,2,2);
imagesc(M1); axis image; axis off;
title('Target texture');
subplot(2,2,3);
imagesc(M1_rec); axis image; axis off;
title('Reconstructed texture');
subplot(2,2,4);
imagesc(E); axis image; axis off;
title('Coordinate Mapping');

%%% perform texture analysis %%%
options.dim = 6;
Jmin = 2;
D = perform_texture_transform(M0,UV,Jmin,+1,options);

return;

%%% save image
rep = ['results/texture_analysis/' name0 '_' name1 '/'];
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
warning on;


[LocP,LocPsi,Segm] = perform_local_pca(Xs, options);