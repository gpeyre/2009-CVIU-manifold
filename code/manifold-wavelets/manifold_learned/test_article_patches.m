% extract and save some local patches from a given texture

rep = 'images/';
name = 'mures';
name = 'dunes';

M = load_image([rep name]);
n = min(128, size(M,1));
M = M(end/2-n/2+1:end/2+n/2,end/2-n/2+1:end/2+n/2,:);
w = 7;
M = rescale(M);

rep = 'results/patches_samples/';
if ~exist(rep)
    mkdir(rep);
end

% save original image
warning off;
imwrite( M, [rep name '_original.png'], 'png' );

p = 1;
k = 0;
while p==1
    k = k+1;
    clf;
    imagesc(M); axis image; axis off;
    [y,x,p] = ginput(1);
    x = round(x);
    y = round(y);
    Mi = M(x-w:x+w,y-w:y+w,:);
    % save image
    imwrite( Mi, [rep name num2str(k) '.png'], 'png' );
end
warning on;