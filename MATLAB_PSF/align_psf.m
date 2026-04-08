function b = align_psf (a, re)
%
%
% b                4D matrix of aligned PSFs
%
% re                 (double)  vector ([x y z]) of gaussian fit radii to
%                              determine subpixel maximum positions.
%                              default is [1 1 3]

if nargin < 1
    error('need input volumes!');
elseif nargin < 2
    re = [1 1 3];
end
% if numel(size(a,3)) == 3
%     sa=size(a);
%     m = zeros(sa(4),3);
%     b = zeros(sa);
%     mid = floor(sa(1:3)/2)+1;
%     [x,y,z]= meshgrid (1:sa(2),1:sa(1),1:sa(3));
%     wbar = waitbar(0,'please wait...');
%     
%     for n=1:sa(4)
%         waitbar (n/(sa(4)),wbar,['aligning n=' num2str(n)]);
%         m(n,:)      = findmax_subpixel(a(:,:,:,n),re) - mid;
%         b(:,:,:,n)  = interp3(x,y,z,a(:,:,:,n),x+m(n,2),y+m(n,1),z+m(n,3),'linear',0);
%     end
%     close(wbar);
% else
%     b = a;
% end
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 500;
b=a;
for i = 1:size(a,4)
b(:,:,:,i) = imregister(a(:,:,:,i),a(:,:,:,1), 'translation', optimizer, metric, 'PyramidLevels', 2);
end