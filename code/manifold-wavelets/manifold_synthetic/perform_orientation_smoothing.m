function Param = perform_orientation_smoothing(Param, options)

% Param = perform_orientation_smoothing(Param, options);

options.null = 0;

Delta = Param.Delta;
Theta = Param.Theta;
n = size(Theta,1);

if isfield(options, 'iter')
    iter = options.iter;
else
    iter = 100;
end
if isfield(options, 'lambda')
    lambda = options.lambda;
else
    lambda = 1;
end

% perform iterative smoothing max(w(:))/2
w = abs(Delta); 
% I = find(w>0); w(I) = max(w(:));
mw = min( max(w(:)), 10 );
w = 1 - min(w/mw,1);
% w = w.^2;
T = Theta;
Tu0 = exp(1i*T);
Tu = Tu0;
h = waitbar(0,'Performing smoothing ...');
    
for i=1:iter
    waitbar(i/iter,h);
    % smooth
    Tu1 = 0*Tu;
    sel1 = mod(-1:n-2,n)+1;
    sel2 = mod(1:n,n)+1;
    Tu1 = (Tu(sel1,:) + Tu(sel2,:) + Tu(:,sel1) + Tu(:,sel2)) / 4;
    % set
    Tu = lambda*Tu1 + (1-lambda)*Tu;
    % assign
    Tu = w.*Tu0 + (1-w).*Tu;
end
close(h);
Param.Theta = mod(angle(Tu),2*pi);
