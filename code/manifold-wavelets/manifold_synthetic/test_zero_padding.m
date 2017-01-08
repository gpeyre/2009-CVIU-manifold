% test for interpolation of edges

n = 512;

Jmin = 2;

w = 11;
options.w = w;
options.n_delta = 11;
options.n_theta = 36;
% sigma = 0.15;
options.manifold_type = 'lines';
options.manifold_type = 'edges';
rep = 'results/';


% wavelet coefficients
ThetaW = rand(2^Jmin)*2*pi;
DeltaW = zeros(2^Jmin);

Jmax = log2(n);

for J=Jmin:Jmax
    ParamW.Theta = zeros(2^J);
    ParamW.Delta = zeros(2^J);
    ParamW.Theta(1:2^Jmin,1:2^Jmin) = ThetaW;
    ParamW.Delta(1:2^Jmin,1:2^Jmin) = DeltaW;
    % synthesis
    % options.sigma = sigma; 
    Param = perform_orientation_transform(ParamW,Jmin,-1,options);
    if J<Jmax-2
        % packing
        w0 = n/2^J + 1;
        options.sigma = sigma*w/w0;
        M = -display_packed_edges(Param,n, options);
    else
        options.sigma = sigma;
        M = perform_edge_reconstruction(Param, options);
    end
    % save image
    warning off;
    imwrite(rescale(M), [rep 'zero_padding_' num2str(J) '.png']);
    warning on;
end