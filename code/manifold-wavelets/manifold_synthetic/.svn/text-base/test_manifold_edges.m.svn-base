% test for edge manifold images

sigma = 0.05*w;
w = 25;
rep = 'results/';

theta_list = [0 0.2 0.4];
delta_list = [-0.3 0 0.3]*w;


x = linspace( -(w-1)/2,(w-1)/2, w );
[Y,X] = meshgrid(x,x);

nt = length(theta_list);
nd = length(delta_list);

k = 0;
for i=1:nt
    for j=1:nd
        k = k+1;
        t = theta_list(i);
        d = delta_list(j);
        y = cos(t)*X + sin(t)*Y - d;
        A = tanh(y/sigma);
        subplot(nd,nt,k);
        imagesc(A);
        axis image; axis off;
        % save image
        warning off;
        imwrite( rescale(A), [rep 'edge_t' num2str(i) '_d' num2str(j) '.png'], 'png' );
    end
end
colormap gray(256);