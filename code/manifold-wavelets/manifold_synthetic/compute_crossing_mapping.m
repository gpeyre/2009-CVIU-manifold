function [delta,theta,match] = compute_crossing_mapping(deltaA,thetaA, deltaB,thetaB, dist, options)

%   [delta,theta,match] = compute_crossing_mapping(deltaA,thetaA,deltaB,thetaB, dist, options);

use_fastc_code = 1;

if ~use_fastc_code
    if length(deltaA(:))>1
        delta = deltaA;
        theta = thetaA;
        match = thetaA;
        for i=1:length(deltaA(:))
            [delta(i),theta(i),match(i)] = compute_edge_mapping(deltaA(i),thetaA(i), deltaB(i),thetaB(i), dist, options);
        end
        return;
    end
end

[n1,n2,s] = size(deltaA);
n = n1*n2; 

% first point
delta1A = deltaA(:,:,1); delta1A = delta1A(:);
theta1A = thetaA(:,:,1); theta1A = mod(theta1A(:),pi);
delta2A = deltaA(:,:,2); delta2A = delta2A(:);
theta2A = thetaA(:,:,2); theta2A = mod(theta2A(:),pi);
% second point
delta1B = deltaB(:,:,1); delta1B = delta1B(:);
theta1B = thetaB(:,:,1); theta1B = mod(theta1B(:),pi);
delta2B = deltaB(:,:,2); delta2B = delta2B(:);
theta2B = thetaB(:,:,2); theta2B = mod(theta2B(:),pi);

if isfield(options, 'w')
    w = options.w;
else
    w = 15;
    options.w = w;
end
w = 15;
options.w = w;

if isfield(options, 'n_delta')
    n_delta = options.n_delta;
else
    n_delta = 13;
    options.n_delta = n_delta;
end
if isfield(options, 'n_theta')
    n_theta = options.n_theta;
else
    n_theta = 20;
    options.n_theta = n_theta;
end

if isfield(options, 'sigma')
    sigma = options.sigma;
else
    sigma = 0.1;
    options.sigma = sigma;
end

dmax = sqrt(2)*w/2+2*sigma*w;
eta = 2;
t = linspace(0,pi-pi/n_theta,n_theta); t = t(:);
d = linspace(-1,1,n_delta);
d = sign(d).*abs(d).^eta; 
d = d(:)*dmax;


%%% retrieve the mappings %%%
global ThetaMapping;
global DeltaMapping;
if isempty(ThetaMapping) || isempty(DeltaMapping)
    if exist('MappingsCrossing.mat')
        load MappingsCrossing;
    end
end
if size(ThetaMapping,1)~=n_delta || size(ThetaMapping,2)~=n_theta
    compute_crossing_mappings(options);
end

% rescale to [0,1]
delta1A = delta1A/dist;
delta1B = delta1B/dist;
delta2A = delta2A/dist;
delta2B = delta2B/dist;

d = d/w;

if ~use_fastc_code
    
    [eps1,idelta1] = min( abs(deltaA-d) );
    [eps2,idelta2] = min( abs(deltaB-d) );

    % [tmp,itheta1] = min( min(abs( [thetaA-t,thetaA-t+2*pi,thetaA-t-2*pi] ), [], 2) );
    % [tmp,itheta2] = min( min(abs( [thetaB-t,thetaB-t+2*pi,thetaB-t-2*pi]
    % ), [], 2) );
    [tmp,itheta1] = min( abs( thetaA-t) );
    [tmp,itheta2] = min( abs( thetaB-t) );

    delta = DeltaMapping(idelta1,itheta1, idelta2,itheta2);
    theta = ThetaMapping(idelta1,itheta1, idelta2,itheta2);
    match = MatchMapping(idelta1,itheta1, idelta2,itheta2);

    % rescale result
    delta = delta*dist;

    if abs(delta)>=dmax*dist/w
        % for un-directional datasets
        theta = compute_circle_mean(thetaA,thetaB);
        delta = sign(delta)*dmax;
        delta = max(abs([deltaA*dist,deltaB*dist,delta]))*sign(delta);
        delta = mean(abs([deltaA*dist,deltaB*dist]))*sign(delta);
        % delta = 6666*sign(delta);
    end

else
   
    d = repmat(d(:)', [n 1]);
    t = repmat(t(:)', [n 1]);
    
    
    [eps1,idelta1A] = min( abs(repmat(delta1A,[1 n_delta])-d), [], 2 );
    [eps2,idelta1B] = min( abs(repmat(delta1B,[1 n_delta])-d), [], 2 );
    [eps1,idelta2A] = min( abs(repmat(delta2A,[1 n_delta])-d), [], 2 );
    [eps2,idelta2B] = min( abs(repmat(delta2B,[1 n_delta])-d), [], 2 );
    
    [eps1,itheta1A] = min( abs(repmat(theta1A,[1 n_theta])-t), [], 2 );
    [eps1,itheta1B] = min( abs(repmat(theta1B,[1 n_theta])-t), [], 2 );
    [eps1,itheta2A] = min( abs(repmat(theta2A,[1 n_theta])-t), [], 2 );
    [eps1,itheta2B] = min( abs(repmat(theta2B,[1 n_theta])-t), [], 2 );
    
    
    DeltaMapping1 = DeltaMapping(:,:,:,:, :,:,:,:, 1);
    DeltaMapping2 = DeltaMapping(:,:,:,:, :,:,:,:, 2);
    ThetaMapping1 = ThetaMapping(:,:,:,:, :,:,:,:, 1);
    ThetaMapping2 = ThetaMapping(:,:,:,:, :,:,:,:, 2);
    
    J = sub2ind(size(DeltaMapping1), idelta1A,itheta1A,idelta2A,itheta2A, idelta1B,itheta1B,idelta2B,itheta2B );
    
    delta = zeros(n1,n2,2); theta = delta;
    delta(:,:,1) = reshape( DeltaMapping1(J), n1, n2);
    delta(:,:,2) = reshape( DeltaMapping2(J), n1, n2);
    theta(:,:,1) = reshape( ThetaMapping1(J), n1, n2);
    theta(:,:,2) = reshape( ThetaMapping2(J), n1, n2);
    match = [];

    % rescale result
    delta = delta*dist;

    I = find( abs(delta)>=dmax*dist/w );
    
    % for un-directional datasets
    theta(I) = compute_circle_mean(thetaA(I),thetaB(I));
    delta(I) = sign(delta(I))*dmax;
    delta(I) = max(abs([deltaA(I)*dist,deltaB(I)*dist,delta(I)]),[],2) .* sign(delta(I));
    delta(I) = sum( abs([deltaA(I)*dist,deltaB(I)*dist]), 2 )/2 .* sign(delta(I));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compute_crossing_mappings(options)

if ~isfield(options, 'n_delta') || ~isfield(options, 'n_theta') || ~isfield(options, 'sigma') || ~isfield(options, 'w')
    error('You must specify options: n_delta, n_theta, sigma, w');
end
n_delta = options.n_delta;
n_theta = options.n_theta;
w = options.w;
w1 = (w-1)/2;
options.rescale = 0;

fprintf('--> Performing edge mapping ...');

[H,delta1_list,theta1_list,delta2_list,theta2_list] = compute_crossing_patches(w,options);
delta1_list = delta1_list/w;
delta2_list = delta2_list/w;

global ThetaMapping;
global DeltaMapping;


[A,D] = compute_association_mapping(H);
ThetaMapping = zeros(n_delta,n_theta,n_delta,n_theta, n_delta,n_theta,n_delta,n_theta, 2);
DeltaMapping = zeros(n_delta,n_theta,n_delta,n_theta, n_delta,n_theta,n_delta,n_theta, 2);
ThetaMapping(:,:,:,:, :,:,:,:, 1) = reshape( theta1_list(A), n_delta,n_theta,n_delta,n_theta, n_delta,n_theta,n_delta,n_theta);
ThetaMapping(:,:,:,:, :,:,:,:, 2) = reshape( theta2_list(A), n_delta,n_theta,n_delta,n_theta, n_delta,n_theta,n_delta,n_theta);
DeltaMapping(:,:,:,:, :,:,:,:, 1) = reshape( delta1_list(A), n_delta,n_theta,n_delta,n_theta, n_delta,n_theta,n_delta,n_theta);
DeltaMapping(:,:,:,:, :,:,:,:, 2) = reshape( delta2_list(A), n_delta,n_theta,n_delta,n_theta, n_delta,n_theta,n_delta,n_theta);

fprintf('done.\n');

save MappingsCrossing *Mapping*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       OLD CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

ThetaMapping = zeros(n_delta,n_theta,n_delta,n_theta);
DeltaMapping = zeros(n_delta,n_theta,n_delta,n_theta);
MatchMapping = zeros(n_delta,n_theta,n_delta,n_theta);

P = zeros(w);

h = waitbar(0,'Computing edge mapping ...');
num = 0;
for d1 = 1:n_delta
for t1 = 1:n_theta
for d2 = 1:n_delta
for t2 = 1:n_theta
    num = num+1;
    waitbar(num/(n_delta*n_delta*n_theta*n_theta),h);
    P1 = H(:,:,d1,t1);
    P2 = H(:,:,d2,t2);
    % compose patch 
    P(1:w1,:) = P1(w1:end,:);
    P(w1:end,:) = P2(1:w1,:);
    P(w1,:) = (P1(end,:)+P2(1,:))/2;
    % perform search
    D = compute_distance_to_points(Hs,P(:));
    [d,i] = min(D);
    DeltaMapping(d1,t1,d2,t2) = delta_list(i);
    ThetaMapping(d1,t1,d2,t2) = theta_list(i);
    MatchMapping(d1,t1,d2,t2) = sqrt(d/w^2);
end
end
end
end
close(h);