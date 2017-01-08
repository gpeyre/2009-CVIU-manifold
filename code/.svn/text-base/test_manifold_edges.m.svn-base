% test for manifold of step edges

n = 200;


rep = 'results/manifold-edges/';
if not(exist(rep))
    mkdir(rep);
end

options.alpha = 5;
M = load_image('fnoise', n, options);

A = cos(rescale(M)*4*pi);
M = A>0;
M = perform_blurring(M,2);

options.bound = 'per';
options.order = 2;
Mh = perform_blurring(M,10);
G0 = grad(Mh, options);
G0 = perform_vf_normalization(G0);

imageplot(M);

niter = 200;
G = G0;
I = find( Mh>.3 & Mh<.7 );
I = [I; I+n^2];
slist = linspace(10,1,niter);
for i=1:niter
    progressbar(i,niter);
    G = perform_blurring(G,slist(i));
    G = perform_vf_normalization(G);
    G(I) = G0(I);
    imageplot(G); drawnow;
end
G = perform_vf_normalization(G);

% angle function
Theta = atan2(G(:,:,2), G(:,:,1));
% distance function
D = perform_redistancing(A, options);

D = D/max(abs(D(:)));
P = cat(1, cos(Theta(:))', sin(Theta(:))', D(:)');

% cut with a straight line

if 0
    x = 1:.5:n; y = {};
    y{1} = (x-n/2)*.3 + n/2-n/4;
    y{2} = (x-n/2)*(-.5) + n/2 + n/4; 
else
    x = {}; y = {};
    for k=1:2
        clf;
        imageplot(G);
        [a,b] = ginput(2);
        t = linspace(0,1,2*n);
        x{k} = a(1)*t+a(2)*(1-t);
        y{k} = b(1)*t+b(2)*(1-t);
    end            
end
c = {'g' 'r' 'b'};

clf;
hold on;
imageplot(M);
for i=1:length(y)
    h = plot(x{i},y{i}, c{i});
    set(h, 'LineWidth', 2);    
end
saveas(gcf, [rep 'cartoon-curves.png'], 'png');

sel = randperm(n^2);
sel = sel(1:min(100^2,n^2));

clf;
hold on;
% draw cylinder
[T,R] = meshgrid(linspace(0,2*pi,100),linspace(-1,1,100));
surf(R,cos(T),sin(T));
shading interp; lighting phong;
colormap gray(256);
% draw points
h = plot3(P(3,sel), P(1,sel), P(2,sel),  'k.');
set(h, 'MarkerSize', 8);
Ind = reshape(1:n^2,n,n);
for i=1:length(y)
    x1 = x{i}; y1 = y{i};
    %    x1 = x1(y1>=1 & y1<=n);
    %    y1 = y1(y1>=1 & y1<=n);
    [Y,X] = meshgrid(1:n,1:n);
    theta = interp2(Y,X,Theta',y1,x1);
    d = interp2(Y,X,D',y1,x1);

    if 1
        s = theta(3:end)-theta(1:end-2);
        s = abs(diff(cos(s)))+abs(diff(sin(s)));
        I = find(s>1)+1; I(end+1) = length(theta);
        prev = 1;
        for k=1:length(I)
            sel = prev:(I(k)-1); prev = I(k)+1;
            h = plot3(d(sel), cos(theta(sel)), sin(theta(sel)), [c{i}]);
            set(h, 'LineWidth', 4);
        end
    else
        h = plot3(d, cos(theta), sin(theta), [c{i}]);
        set(h, 'LineWidth', 4);
    end
end
view(-60,30);
camproj('perspective');
camlight;


axis equal; axis off;
saveas(gcf, [rep 'cartoon-manifold.png'], 'png');

return;

M = double(M>median(M(:)));