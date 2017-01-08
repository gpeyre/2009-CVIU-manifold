function plot_orientation(Param,M, options)


options.null = 0;
if isfield(options, 'do_simplification')
    do_simplification = options.do_simplification;
else
    do_simplification = 0;
end


lw = 1.5;

Theta = Param.Theta;
Delta = Param.Delta;
n = size(Theta,1);


if isfield(options, 'w')
    w = options.w;
else
    w = 11;
end
if isfield(options, 'sigma')
    sigma = options.sigma;
else
    sigma = 0.1;
end

dmax = sqrt(2)*w/2+2*sigma*w;
I = find(abs(Delta)>=dmax);
Delta(I) = Inf * sign(Delta(I));

% remove isolated false orientation
D = abs(Delta);
if do_simplification
    A1 = D(2:end,:);        A1(end+1,:) = Inf;
    A2 = D([1 1:end-1],:);  A2(1,:) = Inf;
    A3 = D(:,2:end);        A3(:,end+1) = Inf;
    A4 = D(:,[1 1:end-1]);  A4(:,1) = Inf;
    I = find( A1==Inf & A2==Inf & A3==Inf & A4==Inf );
    Delta(I) = Inf; D = abs(Delta);
end

vf = zeros(n,n,2);
vf(:,:,1) = -sin(Theta);
vf(:,:,2) = cos(Theta);

if do_simplification
I = find(Delta==Inf | Delta==-Inf);
s = ones(n); s(I) = 0; 
vf(:,:,1) = vf(:,:,1).*s;
vf(:,:,2) = vf(:,:,2).*s;
end

if 0
plot_vf(vf,M');
end

% clf;
x = 1:n;
[Y,X] = meshgrid(x,x);
if do_simplification
    I = find(D~=Inf);
else
    I = (1:n^2)';
end
pos = [X(I), Y(I)];
vect = [-sin(Theta(I)), cos(Theta(I))];
mu = 0.5;
pos1 = pos + vect*mu;
pos2 = pos - vect*mu;

px = [pos1(:,1) pos2(:,1)]';
py = [pos1(:,2) pos2(:,2)]';

hold on;
imagesc(x,x,M');
h = plot( px,py, 'b' );
axis image; axis off;
hold off;
set(h, 'LineWidth', lw);
colormap gray(256);