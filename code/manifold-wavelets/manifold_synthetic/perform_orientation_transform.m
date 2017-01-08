function Param = perform_orientation_transform(Param,Jmin,dir,options)

% perform_orientation_transform - peform non linear wavelet transform
%
% for forward transform
%   ParamW = perform_orientation_transform(Param,Jmin,+1,options);
% or for backward transform
%   Param = perform_orientation_transform(ParamW,Jmin,+1,options);
%
%   Copyright (c) 2006 Gabriel Peyré

n = size(Param.Theta,1);
Jmax = log2(n)-1;

if isfield(options, 'rescale')
    rescale = options.rescale;
else
    rescale = 1;
end

if dir==1
   %%% Forward transform %%%
   Theta = Param.Theta;
   Delta = Param.Delta;
   selx = 1:n;
   sely = 1:n;
   % perform multiscale transform
   for j=Jmax:-1:Jmin
        transfo_axis = [-1,1] * (2*(mod(j,2)==0)-1);
        for s=2:-1:1
            if transfo_axis(s)>0
                % X axis transform
                selxc = 1:selx(end)/2; selyc = sely;
                selxd = selx(end)/2+1:selx(end); selyd = sely;
            else
                % Y axis transform
                selyc = 1:sely(end)/2; selxc = selx;
                selyd = sely(end)/2+1:sely(end); selxd = selx;
            end 
            [Theta(selxc,selyc),Theta(selxd,selyd),Delta(selxc,selyc),Delta(selxd,selyd)] = ...
                    compute_transform_fwd( Theta(selx,sely), Delta(selx,sely), transfo_axis(s), n/2^j, options );
            if rescale
                % rescale
                a = 2^((Jmax-j)+(2-s)/2);
                Theta(selxd,selyd) = Theta(selxd,selyd) * a;
                Delta(selxd,selyd) = Delta(selxd,selyd) * a;
            end
            selx = selxc; sely = selyc;
        end
   end
   Param.Theta = Theta; 
   Param.Delta = Delta;
else
   %%% Backward transform %%%
   ThetaW = Param.Theta; 
   DeltaW = Param.Delta;
   
   selx = 1:2^Jmin;
   sely = 1:2^Jmin;
   
   % perform multiscale transform
   for j=Jmin:Jmax
        transfo_axis = [-1,1] * (2*(mod(j,2)==0)-1);
        for s=1:2
            if transfo_axis(s)>0
                % X axis transform
                selx1 = selx(end)+1:2*selx(end); sely1 = sely;
                selxu = [selx,selx1]; selyu = sely;
            else
                % Y axis transform
                sely1 = sely(end)+1:2*sely(end); selx1 = selx;
                selyu = [sely,sely1]; selxu = selx;
            end 
            if rescale
                a = 2^((Jmax-j)+(2-s)/2);
                ThetaW(selx1,sely1) = ThetaW(selx1,sely1) / a;
                DeltaW(selx1,sely1) = DeltaW(selx1,sely1) / a;
            end
            [ThetaW(selxu,selyu),DeltaW(selxu,selyu)] = compute_transform_bwd( ThetaW(selx,sely), ThetaW(selx1,sely1), DeltaW(selx,sely), DeltaW(selx1,sely1), transfo_axis(s), n/2^j, options );
            selx = selxu; sely = selyu;
        end
   end
   Param.Theta = ThetaW;
   Param.Delta = DeltaW;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ThetaC1,DeltaC1] = compute_transform_bwd( ThetaC, ThetaD, DeltaC, DeltaD, t, d, options )

if isfield(options, 'do_update')
    do_update = options.do_update;
else
    do_update = 1;
end

if t<0
    ThetaC = pi/2 - ThetaC';
    ThetaD = -ThetaD';
    DeltaC = DeltaC';
    DeltaD = DeltaD';
end

% d=spacing between 2 coarse scale coefficients
[n1,n2] = size(ThetaC);

sel1 = 1:n1;        
sel2 = [2:n1,1];  

if do_update
    ThetaC(sel1,:) = ThetaC(sel1,:) - 1/4 * ThetaD;
    ThetaC(sel2,:) = ThetaC(sel2,:) - 1/4 * ThetaD;
    DeltaC(sel1,:) = DeltaC(sel1,:) - 1/4 * DeltaD;
    DeltaC(sel2,:) = DeltaC(sel2,:) - 1/4 * DeltaD;
end

A1 = DeltaC(sel1,:);
A2 = DeltaC(sel2,:);
B1 = ThetaC(sel1,:);
B2 = ThetaC(sel2,:);

% [A,B] = compute_edge_mean(A1,A2,B1,B2,d);
[A,B,match] = compute_edge_mapping(A1,B1,A2,B2,d, options);

%%% perform wavelet mixing %%%
% recover Delta
DeltaC1(1:2:2*n1,:) = DeltaC;
DeltaC1(2:2:2*n1,:) = A + DeltaD;
% recover Theta
ThetaC1(1:2:2*n1,:) = ThetaC;
ThetaC1(2:2:2*n1,:) = B + ThetaD;    
ThetaC1 = mod(ThetaC1,2*pi);

if t<0
    ThetaC1 = pi/2 - ThetaC1';
    DeltaC1 = DeltaC1';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ThetaC,ThetaD,DeltaC,DeltaD] = compute_transform_fwd( Theta, Delta, t, d, options )

if t<0
    Theta = pi/2 - Theta';
    Delta = Delta';
end

if isfield(options, 'do_update')
    do_update = options.do_update;
else
    do_update = 1;
end

Aeven = Delta(1:2:end,:);
Aodd = Delta(2:2:end,:);
Beven = Theta(1:2:end,:);
Bodd = Theta(2:2:end,:);

n1 = size(Aeven,1);
sel1 = 1:n1;        
sel2 = [2:n1,1];
A1 = Aeven(sel1,:);
A2 = Aeven(sel2,:);
B1 = Beven(sel1,:);
B2 = Beven(sel2,:);

% [A,B] = compute_edge_mean(A1,A2,B1,B2,d);
[A,B,match] = compute_edge_mapping(A1,B1,A2,B2,d, options);
B = mod(B,2*pi);

% mean
ThetaC = Beven;
DeltaC = Aeven;
% predict step
ThetaD = mod(Bodd - B,2*pi); 
I = find(ThetaD>pi); ThetaD(I) = ThetaD(I)-2*pi;
DeltaD = Aodd - A;
% update step
if do_update
    ThetaC(sel1,:) = ThetaC(sel1,:) + 1/4 * ThetaD;
    ThetaC(sel2,:) = ThetaC(sel2,:) + 1/4 * ThetaD;
    DeltaC(sel1,:) = DeltaC(sel1,:) + 1/4 * DeltaD;
    DeltaC(sel2,:) = DeltaC(sel2,:) + 1/4 * DeltaD;
end

if t<0
    ThetaC = pi/2 - ThetaC';
    ThetaD = - ThetaD';
    DeltaC = DeltaC';
    DeltaD = DeltaD';
end
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B] = compute_edge_mean(A1,A2,B1,B2,d)

% translate left->right for update
A1trans = A1 - d/2 * cos(B1);
A2trans = A2 + d/2 * cos(B2);
% for angle, no need to translate

%%% compute max interpolation %%%
d1 = abs(A1trans); d2 = abs(A2trans);
if 0
d1 = (cos(B1).*A1 - d/2).^2 + (sin(B1) .* A1).^2;
d2 = (cos(B2).*A2 + d/2).^2 + (sin(B2) .* A2).^2;
end
Bmax = B1;
Amax = A1trans;
I = find( d2<d1 );
Bmax(I) = B2(I);
Amax(I) = A2trans(I);


%%% compute interpolation weight between max/mean %%%
% w=0:orthogonality (take max), 
% w=1:parrallel (take mean)
delta = cos(B1).*cos(B2) + sin(B1).*sin(B2);
w = (delta+1)/2;
% shape this weight
w = reshape_weight(w,1);

%%% compute mean interpolation %%%
% weighted circle mean
if 1
    Bmean = compute_circle_mean(B1,B2);
    % classical mean for 
    Amean = (A1trans + A2trans)/2;
    [Amean,Bmean] = compute_edge_mapping(A1,B1,A2,B2,d);
else
    [Amean,Bmean] = compute_bezier_fit(A1,A2,B1,B2,d,w);
end

%%% compute final mixing %%%
B = compute_circle_mean(Bmean,Bmax,w,1-w);
A = w.*Amean + (1-w).*Amax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B] = compute_bezier_fit(A1,A2,B1,B2,d,w)

if nargin<6
    w = A1*0+1;
end

N = max(4*d, 128);

tension = 1/3;
B = A1*0; A = A1*0;
% use Bezier interpolation to extrapolate
for k=1:prod(size(A1))
    if w(k)>0.05
        
    theta1 = B1(k); r1 = A1(k);
    theta2 = B2(k); r2 = A2(k);
    % compute position of the extremal points
    X1 = [cos(theta1),sin(theta1)] * r1;
    X2 = [d,0] + [cos(theta2),sin(theta2)] * r2;
    s = norme(X1-X2);
    % compute other control points
    X3 = X1 + [sin(theta1),-cos(theta1)]*s*tension;
    X4 = X2 - [sin(theta2),-cos(theta2)]*s*tension;
    d0 = norme(X3-X4);
    X3a = X1 - [sin(theta1),-cos(theta1)]*s*tension;
    X4a = X2 + [sin(theta2),-cos(theta2)]*s*tension;
    da = norme(X3a-X4a);
    swap0 = 1;
    if da<d0
        X3 = X3a; X4 = X4a;
        swap0 = -1;
    end
    % compute the bezier curve
    C = [X1;X3;X4;X2]; % control points
    v = compute_bezier_curve(C,N);
    % find closed point to center [d/2,0]
    dist  = (v(1,:)-d/2).^2 + v(2,:).^2;
    [dist0,I] = min(dist); X = v(:,I); dist0 = sqrt(dist0);
    % assign position and angle
    delta1 = v(:,min(I+1,N)) - v(:,max(I-1,1));
    Btgt = atan2(delta1(2),delta1(1)); % tangent angle
    B(k) = mod(Btgt + pi/2, 2*pi);
    A(k) = dist0;
    % see if we need to swap the directions
    delta2 = [X(1)-d/2, X(2)];
    swap1 = sign( delta1(1)*delta2(2) - delta1(2)*delta2(1) ); 
    % performs the swap
    B(k) = mod( B(k) + (swap0<0)*pi, 2*pi );
    A(k) = A(k)*swap0*swap1;
    % TODO : check that X1 is not too close from X2
    if norme(X1-X2)<1e-2
        %warning('Pbm');
        B(k) = mod( (B1(k)+B2(k))/2,2*pi ); 
        A(k) = ( A1(k) + A2(k) )/2;
    end
    
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = reshape_weight(x,iter)

if nargin<2
    iter = 1;
end
y = 1-(1+cos(pi*x))/2;
for i=1:iter
    y = 1-(1+cos(pi*y))/2;
end