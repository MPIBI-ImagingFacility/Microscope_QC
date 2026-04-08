function cog = centerofgravity(d,thr)
%Determine center of mass of a matrix
% cog = centerofgravity(d,thr)

if nargin == 2
    d(d<thr) = 0;
end
d=double(d);



for dim = 1:numel(size(d))
    dd = d;
    ind = 1:numel(size(d));
    ind(dim) = [];
    dd = permute(dd,[dim ind]);
    dd = reshape(dd,size(d,dim),[]);
    dsum=nansum(dd,2);
    distance = 1:size(d,dim);
    cog(dim) = nansum(dsum(:).*distance(:)) / nansum(dsum);
end