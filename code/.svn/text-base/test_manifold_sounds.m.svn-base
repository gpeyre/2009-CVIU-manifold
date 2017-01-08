% test for spectrogram

name = 'bird';
n = 1024*8;



rep = ['results/manifold/' name '/' ];
if not(exist(rep))
    mkdir(rep);
end
repeps = [rep '/eps/' ];
if not(exist(repeps))
    mkdir(repeps);
end

f = load_sound(name);

% resample
n0 = length(n) / 2;
if strcmp(name, 'bird')
    n0 = length(f) / 4;
end
f = interp1(linspace(0,1,length(f)),f,linspace(0,1,n0)); f = f(:);
% cropn
p = 0;
f = rescale( f(p+(1:n)) );

clf;
plot(f); axis tight;

% windows size
w = 128;
q = 32/2; % overlap
% bell function
x = linspace(0,pi,w);
h = sin(x).^2;

% sampling
xs = 1:q:n-w;
[Y,X] = meshgrid(0:w-1, xs);
sel = X + Y;
H = f(sel)' .* repmat(h(:), [1 size(sel,1)]);
% zero padding
padd = 2;
H(end+1:end+w*(padd-1),:) = 0;
% fft
fH = fft(H,[],1);

% keep only low frequencies
fH = fH(10*padd:end/2,:);


% track max values
[A,R] = max(abs(fH), [], 1);
% treshold to avoid problems
I = find(A<0.3);
R(I) = NaN;
A(I) = NaN;
% same
I = find( abs(diff(R))>10 );
R(I) = NaN;
A(I) = NaN;

% display spectrogram
clf;
plot(f); axis tight;
saveas(gcf, [rep name '-sound.png'], 'png');
saveas(gcf, [repeps name '-sound.eps'], 'epsc');
clf; 
hold on;
imagesc(log(abs(fH))); 
h = plot(R);
axis tight; axis off;
set(h, 'LineWidth', 2);
hold off;
saveas(gcf, [rep name '-spectrogram.png'], 'png');
saveas(gcf, [repeps name '-spectrogram.eps'], 'epsc');

% display spectrogram
clf;
subplot(2,1,1);
plot(f); axis tight;
subplot(2,1,2);
imagesc(log(abs(fH))); axis tight;

% display manifold
clf;
hold on;
h = plot( A, R, 'k' );
set(h, 'LineWidth', 2);
h = scatter(A,R,40, 1:length(A),'filled');
hold off;
axis tight; axis square;
saveas(gcf, [rep name '-manifold.png'], 'png');
saveas(gcf, [repeps name '-manifold.eps'], 'epsc');