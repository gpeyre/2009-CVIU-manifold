function Param = perform_edge_abstraction(M1, options)

% Param = perform_edge_abstraction(M1, options)

options.null = 0;
if isfield(options, 'w')
    w = options.w;
else
    w = 11;
end

n = size(M1,1);
s = w^2;

% the window to weight the distance
if isfield(options, 'bell')
    bell = options.bell;
else
    bell = 'constant';
end

switch lower(bell)
    case 'constant'
        B = ones(w);
    case 'sine'
        x = linspace(0,1,w);
        x = (1-cos(2*pi*x))/2;
        B = x'*x;
    otherwise
        error('Unknown bell shape');
end
B = B(:);


M = extract_subsquares(M1,w);
M = reshape(M,[s,size(M,3)]);
options.rescale = 0;
[H,delta_list,theta_list] = compute_edge_patches(w,options);
H = reshape(H,[s,size(H,3)]);

p = size(M,2);  % #basis functions
q = size(H,2);  % #vectors to parallelize

Delta = zeros(n,n);
Theta = zeros(n,n);

dmax = max(abs(delta_list));

h = waitbar(0,'Computing best fit ...');
for k=1:p
    waitbar(k/p,h);
    x = M(:,k);
    E = sum( (repmat(x,[1,q])-H).^2,1 ); % size 1 x q
    [tmp,I] = min(E);
    Delta(k) = delta_list(I);
    Theta(k) = theta_list(I);
    if abs(Delta(k))==dmax 
        Theta(k) = rand*2*pi;
    end       
    if abs(Delta(k))==dmax && 0
        mu = [];
        if k-1>0
            mu = [mu Theta(k-1)];
        end
        if k-n+1>0
            mu = [mu Theta(k-n+1)];
        end
        if k-n>0
            mu = [mu Theta(k-n)];
        end
        if k-n-1>0
            mu = [mu Theta(k-n-1)];
        end
        if ~isempty(mu)
            Theta(k) = mean(mu);
        end
    end
end
close(h);

Param.Delta = Delta;
Param.Theta = Theta;


return;


La = zeros(n,n);
Mu = zeros(n,n);

% <1,1>_B
Bs = sum(B); 

% <M,1>_B
m = sum(M .* repmat(B,[1,p]),1); % mean
% <H,1>_B
B = repmat(B,[1,q]);
ct = sum(H.*B,1); % size 1 x q

h = waitbar(0,'Computing best fit ...');
for k=1:p
    waitbar(k/p,h);
    x = M(:,k);
    % <x_k,phi_n> for all n
    mt = sum( repmat(x,[1,q]).*H.*B,1 ); % size 1 x q
    mu = (mt.*ct-m(k)) ./ (ct.^2-Bs);
    la = (Bs*mt-m(k)*ct) ./ (Bs-ct.^2);
    E = repmat(la,s,1).*H + repmat(mu,s,1) - repmat(x,1,q);
    % |M-la*H+mu|_B^2
    E = sum( B.*(E).^2,1 ); % error
    I = find(la<0); E(I) = Inf;
    [tmp,I] = min(E);
    Delta(k) = delta_list(I);
    Theta(k) = theta_list(I);
    La(k) = la(I);
    Mu(k) = mu(I);
end
close(h);

Param.Delta = Delta;
Param.Theta = Theta;
Param.Lambda = La;
Param.Mu = Mu;