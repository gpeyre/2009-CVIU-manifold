% test for BW image crossing synthesis

n = 128;

Jmin = 2;

options.w = 15;
options.n_delta = 5;
options.n_theta = 8;

% options.sigma = 0.1;

% wavelet coefficients
ThetaW = zeros(n,n,2);
DeltaW = zeros(n,n,2);
ThetaW(1:2^Jmin,1:2^Jmin,1:2) = rand(2^Jmin,2^Jmin,2)*2*pi;

options.null = 0;
ParamW.Theta = ThetaW;
ParamW.Delta = DeltaW;

% synthesis
Param = perform_crossing_transform(ParamW,Jmin,-1,options);

return;

M1 = perform_edge_reconstruction(Param, options);

clf;
imagesc(M1);
colormap gray(256);

return;

clf;
plot_orientation(Param,M1, options);
colormap gray(256);