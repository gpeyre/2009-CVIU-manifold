function [delta,theta,match] = compute_edge_mapping(delta1,theta1, delta2,theta2, dist, options)

%   [delta,theta,match] = compute_edge_mapping(delta1,theta1,delta2,theta2, dist, options);

use_fastc_code = 1;

if ~use_fastc_code
    if length(delta1(:))>1
        delta = delta1;
        theta = theta1;
        match = theta1;
        for i=1:length(delta1(:))
            [delta(i),theta(i),match(i)] = compute_edge_mapping(delta1(i),theta1(i), delta2(i),theta2(i), dist, options);
        end
        return;
    end
end

[n1,n2] = size(delta1);
n = n1*n2; 
delta1 = delta1(:);
delta2 = delta2(:);
theta1 = theta1(:);
theta2 = theta2(:);


theta1 = mod(theta1,2*pi);
theta2 = mod(theta2,2*pi);

if isfield(options, 'w')
    w = options.w;
else
    w = 15;
    options.w = w;
end
w = 27;
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
t = linspace(0,2*pi-2*pi/n_theta,n_theta); t = t(:);
d = linspace(-1,1,n_delta);
d = sign(d).*abs(d).^eta; 
d = d(:)*dmax;


%%% retrieve the mappings %%%
global ThetaMapping;
global DeltaMapping;
global MatchMapping;
if isempty(ThetaMapping) || isempty(DeltaMapping)
    if exist('Mappings.mat')
        load Mappings;
    end
end
if size(ThetaMapping,1)~=n_delta || size(ThetaMapping,2)~=n_theta
    compute_mappings(options);
end

% rescale to [0,1]
delta1 = delta1/dist;
delta2 = delta2/dist;
d = d/w;

if ~use_fastc_code
    
    [eps1,idelta1] = min( abs(delta1-d) );
    [eps2,idelta2] = min( abs(delta2-d) );

    % [tmp,itheta1] = min( min(abs( [theta1-t,theta1-t+2*pi,theta1-t-2*pi] ), [], 2) );
    % [tmp,itheta2] = min( min(abs( [theta2-t,theta2-t+2*pi,theta2-t-2*pi] ), [], 2) );
    [tmp,itheta1] = min( abs( theta1-t) );
    [tmp,itheta2] = min( abs( theta2-t) );

    delta = DeltaMapping(idelta1,itheta1, idelta2,itheta2);
    theta = ThetaMapping(idelta1,itheta1, idelta2,itheta2);
    match = MatchMapping(idelta1,itheta1, idelta2,itheta2);

    % rescale result
    delta = delta*dist;

    if abs(delta)>=dmax*dist/w
        % for un-directional datasets
        theta = compute_circle_mean(theta1,theta2);
        delta = sign(delta)*dmax;
        delta = max(abs([delta1*dist,delta2*dist,delta]))*sign(delta);
        delta = mean(abs([delta1*dist,delta2*dist]))*sign(delta);
        % delta = 6666*sign(delta);
    end

else
   
    d = repmat(d(:)', [n 1]);
    t = repmat(t(:)', [n 1]);
    
    [eps1,idelta1] = min( abs(repmat(delta1,[1 n_delta])-d), [], 2 );
    [eps2,idelta2] = min( abs(repmat(delta2,[1 n_delta])-d), [], 2 );

    [eps1,itheta1] = min( abs(repmat(theta1,[1 n_theta])-t), [], 2 );
    [eps1,itheta2] = min( abs(repmat(theta2,[1 n_theta])-t), [], 2 );
    
    
    J = sub2ind(size(DeltaMapping), idelta1,itheta1, idelta2,itheta2 );

    delta = DeltaMapping(J);
    theta = ThetaMapping(J);
    match = MatchMapping(J);

    % rescale result
    delta = delta*dist;

    I = find( abs(delta)>=dmax*dist/w );
    
    % for un-directional datasets
    theta(I) = compute_circle_mean(theta1(I),theta2(I));
    delta(I) = sign(delta(I))*dmax;
    delta(I) = max(abs([delta1(I)*dist,delta2(I)*dist,delta(I)]),[],2) .* sign(delta(I));
    delta(I) = sum( abs([delta1(I)*dist,delta2(I)*dist]), 2 )/2 .* sign(delta(I));
    
    delta = reshape(delta,n1,n2);
    theta = reshape(theta,n1,n2);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compute_mappings(options)

if ~isfield(options, 'n_delta') || ~isfield(options, 'n_theta') || ~isfield(options, 'sigma') || ~isfield(options, 'w')
    error('You must specify options: n_delta, n_theta, sigma, w');
end
n_delta = options.n_delta;
n_theta = options.n_theta;
w = options.w;
w1 = (w-1)/2;
options.rescale = 0;
p = n_delta*n_theta; %  number of patches

fprintf('--> Performing edge mapping ...');

[H,delta_list,theta_list] = compute_edge_patches(w,options);
delta_list = delta_list/w;

global ThetaMapping;
global DeltaMapping;
global MatchMapping;

[A,D] = compute_association_mapping(H);
ThetaMapping = reshape( theta_list(A), n_delta,n_theta,n_delta,n_theta);
DeltaMapping = reshape( delta_list(A), n_delta,n_theta,n_delta,n_theta);
MatchMapping = reshape( D, n_delta,n_theta,n_delta,n_theta);;

fprintf('done.\n');

save Mappings *Mapping*

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