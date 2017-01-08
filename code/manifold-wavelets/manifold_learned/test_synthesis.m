% test for manifold based synthesis

path(path,'nntools');

if ~exist('name')
    name = 'sea_grid';
    name = 'corral';
    name = 'framboise';
    name = 'mures';
end

rep = 'images/';

p = 64;
M0 = load_image([rep name]);
n0 = 64;
n0 = min(min(size(M0,1),size(M0,2)),n0);
M0 = M0(end/2-n0/2+1:end/2+n0/2, end/2-n0/2+1:end/2+n0/2,:);
M0 = rescale(M0);
options.wmax = 4;

% generate random coarse coefficients
Jmin = 1;
n1 = 64;
D = zeros(n1,n1,2);
Coarse = floor(rand(2^Jmin,2^Jmin,2)*n0)+1;
D(1:2^Jmin,1:2^Jmin,:) = Coarse;


options.null = 0;
[UVfinal,UVint] = perform_texture_transform(M0,D,Jmin,-1,options);

rep = ['./results/texture_synthesis/' name '/' ];
for i=3:length(UVint)
    UV = UVint{i};
    m = size(UV,1);
    % perform reconstruction
    M1 = perform_texture_reconstruction(M0,UV);

    % the mapping
    E = zeros(m,m,3);
    E(:,:,1:2) = rescale(UV);

    % save results
    if ~exist(rep)
        mkdir(rep);
    end
    warning off;
    imwrite( M0, [rep name '_original.png'], 'png' );
    imwrite( M1, [rep name '_synthetized_' num2str(m) '.png'], 'png' );
    imwrite( E, [rep name '_mapping_' num2str(m) '.png'], 'png' );
    warning on;
end


% display results
clf;
subplot(1,2,1);
d = (n1-n0)/2;
imagesc(M0); axis image; axis([1-d n1-d 1-d n1-d]); axis off;
title('Original');
subplot(1,2,2);
imagesc(M1); axis image; axis off;
title('Synthetised');
if size(M0,3)==1
    colormap gray(256);
end
rep = './results/texture_synthesis/';
saveas(gcf, [rep name '_synthesis.png'], 'png');