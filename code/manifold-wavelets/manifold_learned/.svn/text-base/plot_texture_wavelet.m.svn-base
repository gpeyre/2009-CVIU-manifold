function Mout = plot_texture_wavelet(M, Jmin, options)

% plot_texture_wavelet - plot 2D non linear wavelet transform
%
%   Mout = plot_texture_wavelet(M, Jmin, options);
%
%   'MW' is the wavelet transform (in Mallat's ordering, not lifting inplace
%       ordering).
%   'Jmin' is the minimum scale of the transform.
%   
%   Copyright (c) 2006 Gabriel Peyré

if nargin<2
    Jmin = 1;
end

options.null = 0;
if isfield(options, 'style')
    style = options.style;
else
    style = 'real';
end
if isfield(options, 'edge_color')
    edge_color = options.edge_color;
else
    edge_color = 'r';
end
if isfield(options, 'style')
    renormalize = options.renormalize;
else
    renormalize = 0;
end
if isfield(options, 'line_width')
    lw = options.line_width;
else
    lw = 2;
end

M = sqrt(sum(M.^2,3));

n = size(M,1);
Jmax = log2(n)-1;

for j=Jmin:Jmax
    qmin = ~(j==Jmin);
    for q=qmin:2
        [selx,sely] = perform_subband_selection(j,q);
        M1 = M(selx,sely,:);

        if 0 % q>0
            M1 = 0.3 * M1/std(M1(:));
            M1 = rescale( clamp(M1, -1,1) );
        else
            M1 = rescale(M1);
        end
            
        M(selx,sely,:) = M1;
    end
end

hold on;
s = [1/2 n-1/2];
imagesc(s,s,M);
axis image; axis off;
colormap gray(256);
% display axis separation
for j=Jmax:-1:Jmin;
    x = [0 2^(j+1)];
    y = [2^j 2^j];
    if mod(j,2)==1
        h = plot(y,x, edge_color);
    else
        h = plot(x,y, edge_color);
    end
    set(h, 'LineWidth', lw);
    y = [0 2^(j)];
    x = [2^j 2^j];
    if mod(j,2)==1
        h = plot(y,x, 'r');
    else
        h = plot(x,y, 'r');
    end
    set(h, 'LineWidth', lw);
end
% plot overall square
x = [0 0 n n 0];
y = [0 n n 0 0];
h = plot(x,y,edge_color);
set(h, 'LineWidth', lw*1.2);
hold off;
axis ij;

if nargout>1
    Mout = M;
end