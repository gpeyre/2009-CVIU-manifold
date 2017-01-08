% test for shape animation

n = 256;

Jmin = 3;

load Mappings;
options.w = 15;
options.n_delta = size(DeltaMapping,1);
options.n_theta = size(DeltaMapping,2);

% sigma = 0.15;
options.manifold_type = 'lines';
options.manifold_type = 'edges';
rep = 'results/';

% wavelet coefficients
ThetaW = rand(2^Jmin)*2*pi;
DeltaW = zeros(2^Jmin);

Jmax = log2(n);

nbr_frames = 100;
rep = 'results/anim/';
if ~exist(rep)
    mkdir(rep);
end

ThetaOffs = ones(2^Jmin)*2*pi/nbr_frames;
ThetaOffs = ThetaOffs .* sign( rand(2^Jmin)-0.5 );
if 0
mask = zeros(2^Jmin);
mask(2,2) = 1;
ThetaOffs = ThetaOffs .* mask;
end
h = waitbar(0,'Computing animation.');
for i=1:nbr_frames
    waitbar(i/nbr_frames,h);
    % synthesis
    ParamW.Theta = zeros(n);
    ParamW.Delta = zeros(n);
    ParamW.Theta(1:2^Jmin,1:2^Jmin) = ThetaW + ThetaOffs*(i-1);
    ParamW.Delta(1:2^Jmin,1:2^Jmin) = DeltaW;
    % options.sigma = sigma; 
    Param = perform_orientation_transform(ParamW,Jmin,-1,options);
    M = perform_edge_reconstruction(Param, options);
    % save image
    warning off;
    imwrite(rescale(M), [rep 'anim_' num2string_fixeddigit(i, 3) '.png']);
    warning on;
end
close(h);