% test for edges transform

n = 128;
options.w = 21;
options.n_delta = 9;
options.n_theta = 12;
options.sigma = 0.1;
options.bell = 'constant';
options.bell = 'sine'; % to concentrate the error on the center

if 1
    name = 'disk';
    name = 'square';
    name = 'polygons_blurred';
    M = load_image(name, n);
else
    name = 'lena';
    name = 'peppers';
    M = load_image(name,256);
    M = M(end/2-n/2+1:end/2+n/2, end/2-n/2+1:end/2+n/2);
end

M = rescale(M);

Param = perform_edge_abstraction(M,options);

options.sigma = 0.1;
M1 = perform_edge_reconstruction(Param, options);

clf;
subplot(1,2,1);
if 0
end
imagesc(M);
axis image; axis off;
colormap gray(256);

subplot(1,2,2);
imagesc(M1);
title('Reconstructed.');
axis image; axis off;

return;

% test for wavelet transform
ParamW = perform_edge_wavelet_transform(Param, +1, options);

% perform thresholding of the distance map
T = 50; 
ParamWT = ParamW;
% ParamWT.Full = keep_above(ParamWT.Full, T);
ParamWT.Theta = keep_above(ParamWT.Theta, T);
Param1 = perform_edge_wavelet_transform(ParamWT, -1, options);
MT = perform_edge_reconstruction(Param1, options);

clf;
subplot(1,2,1);
imagesc(M1); 
axis image; axis off;
subplot(1,2,2);
imagesc(MT); 
axis image; axis off;
colormap gray(256);
