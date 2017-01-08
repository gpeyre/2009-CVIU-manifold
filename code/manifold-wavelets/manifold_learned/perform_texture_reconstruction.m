function M1 = perform_texture_reconstruction(M0,UV,ntgt)

%   M1 = perform_texture_reconstruction(M0,UV,ntgt);

s = size(M0,3);

U = UV(:,:,1);
V = UV(:,:,2);

n1 = size(U,1);
if nargin<3
    ntgt = n1;
end

I = sub2ind(size(M0),U,V);
if n1>=ntgt
    
    M1 = zeros(n1,n1,s);
    for i=1:s
        Mu = M0(:,:,i);
        M1(:,:,i) = Mu(I);
    end

else
    % perform extrapolation
    w = round( ntgt/(2*n1) );
    ntgt = n1 * 2 * w;
    ww = 2*w+1;
    M1 = zeros(ntgt,ntgt,s);;
    options.sampling = 'uniform';
    X0 = compute_patch_library(M0,w, options);
    k = 0;
    r = ntgt/n1;
    for i=1:n1
        for j=1:n1
            k = k+1;
            selx = (i-1)*r+1:(i-1)*r+ww-1;
            sely = (j-1)*r+1:(j-1)*r+ww-1;
            selx = mod(selx-1,ntgt)+1;
            sely = mod(sely-1,ntgt)+1;
            M1(selx,sely,:) = X0(1:end-1,1:end-1,:,I(k));
        end
    end
end