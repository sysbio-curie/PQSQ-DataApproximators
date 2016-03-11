function [xsorted_infs,indices,intervals_symmetrical] = prepareForSplitPointIntervalsFast(x,intervals)

[xsorted,indices] = sort(x);

intervals_symmetrical = makeSymmetricIntervals(intervals);

xsorted_infs = zeros(size(x,1)+2,1);

xsorted_infs(1) = -Inf;
for i=1:size(x,1)
    xsorted_infs(i+1) = xsorted(i);
end
xsorted_infs(size(xsorted_infs,1)) = Inf;