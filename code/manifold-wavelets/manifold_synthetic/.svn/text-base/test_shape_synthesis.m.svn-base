% test for shape synthesis

n = 128*2;
clear options;
options.sigma = 0.05;
options.w = 9;
options.n_delta = 9;
options.n_theta = 12;


save_images = 1;
rep = 'results/';


Jmin = 1;
Jmax = log2(n)-1;
ParamW.Theta = zeros(n);
ParamW.Delta = zeros(n);
ParamW.Theta(1:2^Jmin,1:2^Jmin) = rand(2^Jmin) * 2*pi;
ParamW.Delta(1:2^Jmin,1:2^Jmin) = 0;

thetamax = 2;
deltamax = 5;
alpha_list = linspace(5,0.01,9);
sigma_min = 0.05;

M = {};
for alpha = alpha_list
    fprintf('.');
    % compute decreasing spectra
    nrb_j = (Jmax-Jmin+1)*2;
    sigma = linspace(1,sigma_min,nrb_j);
	sigma_theta = thetamax * sigma.^alpha;
    sigma_delta = deltamax * sigma.^alpha;
    % set up coefficients
    for j = Jmin:Jmax
        for s = 1:2
            [selx,sely] = perform_subband_selection(j,s); 
            k = (j-Jmin)*2+s;
            ParamW.Theta(selx,sely) = randn(length(selx),length(sely)) * sigma_theta( k );
            ParamW.Delta(selx,sely) = randn(length(selx),length(sely)) * sigma_delta( k );
        end
    end
    Param = perform_orientation_transform(ParamW,Jmin,-1,options);
    M{end+1} = perform_edge_reconstruction(Param, options);
end
fprintf('\n');

%%% display %%%
p = ceil(sqrt(length(M)));
clf;
warning off;
for i=1:length(M)
    subplot(p,p,i);
    imagesc(M{i});
    axis image; axis off;
    % save image
    imwrite(rescale(M{i}), [rep 'shape_synthesis_' num2str(i) '.png'], 'png');
end
colormap gray(256);
warning on;