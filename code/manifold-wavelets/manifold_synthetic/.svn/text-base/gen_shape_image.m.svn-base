function [M,Param] = gen_shape_image(n, options)

ampl = 8;
n = ampl*n;

options.null = 1;
if ~isfield(options, 'alpha')
    options.alpha = 0.8;
end   

name = 'rand';
p = n * 80;
c = load_curve(name, p, options);

c(:,end) = [];

s = 0.1;
c(1,:) = rescale(c(1,:),s,1-s)*n + 1;
c(2,:) = rescale(c(2,:),s,1-s)*n + 1;

c = perform_curve_resampling(c, 0.2, 'dist');
p = size(c,2);

Ix = round(c(1,:));
Iy = round(c(2,:));

A = zeros(n); % image

for i=1:p
    A(Ix(i),Iy(i)) = 1;
end

%%% compute M %%%
M = bwlabel(1-A,4);
M = double(M==M(end/2,end/2));
M = 2*M-1;

%%% compute Delta and Theta %%%
% distance transform
Delta = M*0;
[D,L] = bwdist(M==1);
I = find(D>0);
Delta(I) = D(I)-1/2;
[D,L] = bwdist(M==-1);
I = find(D>0);
Delta(I) = -D(I)+1/2;
% smooth a bit
h = ones(6);
h = conv2(h,h);
h = h/sum(h(:));
Delta1 = perform_convolution(Delta,h);
[gx,gy] = compute_grad(Delta1);
Theta = atan2(gy,gx);


%%% assign %%%
Param.Theta = Theta;
Param.Delta = Delta;

%%% sub sample %%%
s = floor(ampl/2)*2+1;
h = ones(s); h  = h/sum(h(:));
Param.Delta = perform_convolution(Param.Delta, h);
% Param.Theta = angle( perform_convolution( exp(1i*Param.Delta), h) );

% subsample
M = imresize(M,1/ampl,'bicubic');
Param.Delta = imresize(Param.Delta,1/ampl,'bicubic');
Param.Delta = Param.Delta/ampl;
Param.Theta = imresize( exp(1i*Param.Theta),1/ampl,'bicubic');
Param.Theta = angle( Param.Theta );