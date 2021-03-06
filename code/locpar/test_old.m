
return;



% detect the central frequency
HF = fftshift( abs(fft2(H)) );
HF = mean(HF,4);
HF(end/2+1,end/2+1) = 0;
[Y,X] = meshgrid(-w/2:w/2-1,-w/2:w/2-1);
theta = atan2(Y,X);
r = sqrt(X.^2+Y.^2);

plot(r(:), HF(:), '.');

return;

options.order = 2;
G = grad( perform_blurring(M,1.3) ,options);
G = G(:,:,2:-1:1); G(:,:,1) = -G(:,:,1);

w = 16*2;
options.sub = 8;
H = compute_all_patch(M, w, options);
Hg = cat(3, compute_all_patch(G(:,:,1), w, options), ...
            compute_all_patch(G(:,:,2), w, options) );
T = cat(3,Hg(:,:,1,:).^2,Hg(:,:,2,:).^2, prod(Hg,3));
T = squeeze( sum(sum(T)) );
p = size(T,2);
T = reshape(T', [sqrt(p) sqrt(p) 3]);
[e1,e2,l1,l2] = perform_tensor_decomp(T);
U = perform_tensor_mapping(T,+1);
Theta = U(:,:,3);

% number of patches
q = 1000;
sel = randperm(p); sel = sel(1:q);
[Y,X] = meshgrid(1:w,1:w);
X = repmat(X,[1 1 q]);
Y = repmat(Y,[1 1 q]);
ThetaS = repmat( reshape(Theta(sel),[1 1 q]), [w w 1]);
% sampling locations
S = sin(ThetaS).*X - cos(ThetaS).*Y;
S = reshape(S, [w*w q]);
HS = reshape(H(:,:,sel), [w*w q]);

return;


% compute a profile by cliking
w = 40; % length of the profile
b = [1 1];
p = 100;
L = [];
while true;
    clf;
    imagesc(M); axis image; axis off;
    colormap gray(256);
    [y,x,b] = ginput(2);
    if b(1)~=1 || b(2)~=1
        break;
    end
    p1 = [x(1) y(1)]; p2 = [x(2) y(2)];
    d = norm(p1-p2);
    if d>w
        t = linspace(0,w/d,p);        
        v = interp2(M, x(1)*t+x(2)*(1-t), y(1)*t+y(2)*(1-t));
        L = [L, v(:)];
    else
        warning('Too short');
    end
end

return;