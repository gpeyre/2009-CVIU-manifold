% test for direction estimation + smoothing

n = 64;
clear options;
sigma = 0.1;
options.sigma = sigma;
options.rescale = 0;

if 0
    options.alpha = 3;  
    M = load_image('fnoise', n, options);
    M = 2*double( rescale(M,-1,1)>0 ) - 1;
    options.w = 7;
    options.n_delta = 15;
    options.n_theta = 20;
    Param = perform_edge_abstraction(M, options);
else
    options.alpha = 0.9;
    [M,Param] = gen_shape_image(n, options);
    dmax = 6;
    I = find(abs(Param.Delta)>dmax); Param.Theta(I) = rand(length(I),1)*2*pi;
    I = find(Param.Delta>dmax); Param.Delta(I) = dmax;
    I = find(Param.Delta<-dmax); Param.Delta(I) = -dmax;
end

h = compute_gaussian_filter([1 1]*13,8*sigma/n,[n n]);
M = perform_convolution(M,h);

save_images = 1;
rep = 'results/binary_shapes/';
if ~exist(rep)
    mkdir(rep);
end

M1 = perform_edge_reconstruction(Param, options);

options.iter = 1000;
options.lambda = 1;
Param1 = perform_orientation_smoothing(Param, options);

if save_images && n==64
    clf;
    subplot(1,2,1);
    plot_orientation(Param,M, options);
    subplot(1,2,2);
    plot_orientation(Param1,M, options);
    saveas(gcf, [rep 'orientation_smoothing.eps'], 'epsc');
    % save distance field
    warning off;
    imwrite( rescale(Param.Delta), [rep 'distance.png'], 'png' );
    warning on;
    A = Param.Delta<dmax & Param.Delta>-dmax;
    imwrite( 1-rescale(A)*0.6, [rep 'band.png'], 'png' );
end

options.do_update = 1;
options.n_delta = 9;
options.n_theta = 12;
% try some wavelet transform
Jmin = 2;
ParamW = perform_orientation_transform(Param1,Jmin,+1,options);
Param2 = perform_orientation_transform(ParamW,Jmin,-1,options);

if save_images
    warning off;
    imwrite(rescale(M), [rep 'shape.png'], 'png');
    imwrite(rescale(Param1.Delta), [rep 'delta.png'], 'png');
    imwrite(rescale(cos(Param1.Theta)), [rep 'theta_cos.png'], 'png');
    imwrite(rescale(sin(Param1.Theta)), [rep 'theta_sin.png'], 'png');
    warning on;
    plot_orientation_wavelet(ParamW.Theta, Jmin);
    saveas(gcf, [rep 'theta_wav.eps'], 'epsc');
    plot_orientation_wavelet(ParamW.Delta, Jmin);
    saveas(gcf, [rep 'delta_wav.eps'], 'epsc');
end

% test for bijectivity
D = mod(Param1.Theta-Param2.Theta,2*pi);
I = find(D>pi); D(I) = D(I)-2*pi;
err = norm(Param1.Delta-Param2.Delta) + norm(D);
disp(['Error (should be 0) ' num2str(err,2) '.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for linear approximation
Jmax = log2(n)-1;
nbrJ = 9; a = ceil(sqrt(nbrJ));
clf;
for Jmin1 = Jmax+1:-1:max(Jmax+2-nbrJ,2)
    fprintf('.');
    % keep low frequency content
    ParamW1.Theta = zeros(n); 
    ParamW1.Theta(1:2^Jmin1,1:2^Jmin1) = ParamW.Theta(1:2^Jmin1,1:2^Jmin1);
    ParamW1.Delta = zeros(n); 
    ParamW1.Delta(1:2^Jmin1,1:2^Jmin1) = ParamW.Delta(1:2^Jmin1,1:2^Jmin1);
    % reconstruct
    Param2 = perform_orientation_transform(ParamW1,Jmin,-1,options);
    M1 = perform_edge_reconstruction(Param2, options);
    subplot(a,a,Jmax+2-Jmin1);
    imagesc(M1); axis image; axis off;
    title(['J=' num2str(Jmax+1-Jmin1)]);
    if save_images
        warning off;
        imwrite(rescale(M1), [rep 'approx_lin_j' num2str(Jmin1) '.png'], 'png');
        warning on;
    end
end
fprintf('\n');
colormap gray(256);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for NL approximation
Jmax = log2(n)-1;
T_theta = [5 10 15 20 25 30];
T_delta = [0.1 0.2 0.6 0.8 1];
nbr = length(T_theta); a = ceil(sqrt(nbrJ));
clf;
for i=1:nbr
    fprintf('.');
    ParamW1.Theta = keep_above(ParamW.Theta,T_delta(i)); 
    ParamW1.Delta = keep_above(ParamW.Delta,T_theta(i)); 
    % reconstruct
    Param2 = perform_orientation_transform(ParamW1,Jmin,-1,options);
    M1 = perform_edge_reconstruction(Param2, options);
    subplot(a,a,i);
    imagesc(M1); axis image; axis off;
    title(['J=' num2str(Jmax+1-Jmin1)]);
end
fprintf('\n');
colormap gray(256);
