function  [I, V, E] = perform_kmeans_distance(D, k, options)

% number of points
n = size(D,1);

options.null = 0;
if isfield(options, 'niter')
    niter = options.niter;
else
    niter = 20;
end

if isfield(options, 'nrestart')
    nrestart = options.nrestart;
else
    nrestart = 1;
end

if nrestart>1
    options.nrestart = 1;
    emin = Inf;
    for i=1:nrestart
        [I1,V1,E1] = perform_kmeans_distance(D, k, options);
        if E1(end)<emin
            emin = E1(end);
            I = I1; V = V1; E = E1;
        end
    end
    return;
end

if isfield(options, 'V')
    V = options.V;
else
    % random starting points by farthest sampling
    V(1) = floor( rand*(n-1) ) + 1;
    V(1) = pick_farthest(D,V);
    for i=2:k
        V(end+1) = pick_farthest(D,V);
    end
end

% replace distances by square distances
D = D.^2;
E = [];
for i=1:niter
    progressbar(i,niter);
    % compute the partition function
    [tmp,I] = min( D(V,:) );
    % find barycenters
    E(end+1) = 0;
    for j=1:k
        J = find(I==j);
        if not(isempty(J))
            [e,V(j)] = min( sum( D(J,:) ) );
            E(end) = E(end) + e;
        end
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = pick_farthest(D,V)

D = min( D(:,V), [], 2  );
[tmp,V] = max( D );
