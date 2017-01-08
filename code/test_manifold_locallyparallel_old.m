% test for spatially varying kernel defined over steerable wavelets

n = 128*2;

name = 'barb';


rep = ['results/manifold/locallypar/' ];
if not(exist(rep))
    mkdir(rep);
end


M = load_image(name);
M = crop(M,n);
M = rescale(sum(M,3));

w = 17; q = 8;
w = 11; q = 5;
w = 11; q = 2;


warning off;
imwrite(M, [rep name '-original.png'], 'png');
warning off;


% test for local fourier mode
niter = 5;
x = []; y = [];
for i=1:niter
    options.window_type = 'constant';
    options.window_type = 'sin';
    [MF,x(end+1),y(end+1)] = grab_fourier_mode(M,w);
    warning off;
    imwrite(rescale(abs(MF)), [rep name '-locfourier-' num2str(i) '.png'], 'png');
    warning off;
end


clf;
hold on;
imagesc(M); axis image; axis off;
plot(y,x, '.');
hold on;
axis ij;
saveas(gcf, [rep name '-points.png'], 'png');

