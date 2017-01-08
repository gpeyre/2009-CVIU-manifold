% test for BW image synthesis

n = 512;

Jmin = 5;

options.w = 15;
options.n_delta = 23;
options.n_theta = 80;

% options.sigma = 0.1;
options.manifold_type = 'lines';
options.manifold_type = 'edges';

% wavelet coefficients
ThetaW = zeros(n);
DeltaW = zeros(n);
ThetaW(1:2^Jmin,1:2^Jmin) = rand(2^Jmin)*2*pi;
% DeltaW(1:2^Jmin,1:2^Jmin) = randn(2^Jmin)*40;
% ThetaW(1:2^Jmin,1:2^Jmin) = [1,-1;-1,1]*pi/4;
e = 0;
%ThetaW(1:2^Jmin,1:2^Jmin) = [ [pi/2+e;-pi/2+e] [pi/2-e;-pi/2-e] ]+pi/2;
% ThetaW(1:2^Jmin,1:2^Jmin) = [ [e;-e] [pi+e;pi-e] ]+pi/2;

%ThetaW(1:2^Jmin,1:2^Jmin) = [ [e;pi/2-e] [pi/2+e;e] ]+pi/2;
% ThetaW(1:2^Jmin,1:2^Jmin) = ThetaW(1:2^Jmin,1:2^Jmin) + randn(2^Jmin)*0.1;

% DeltaW(1:2^Jmin,1:2^Jmin) = 5;
% DeltaW(1:2^Jmin,1:2^Jmin) = [ [-5;5] [-5;5] ];
options.null = 0;


ParamW.Theta = ThetaW;
ParamW.Delta = DeltaW;

% synthesis
Param = perform_orientation_transform(ParamW,Jmin,-1,options);
M1 = perform_edge_reconstruction(Param, options);

clf;
imagesc(M1);
colormap gray(256);

return;

clf;
plot_orientation(Param,M1, options);
colormap gray(256);