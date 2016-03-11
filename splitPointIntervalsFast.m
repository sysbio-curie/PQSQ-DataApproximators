function inds = splitPointIntervalsFast(xsorted_infs,indices,beta,intervals_symmetrical)

nint = (size(intervals_symmetrical,2))/2;

ints1 = intervals_symmetrical+beta;

h = histc(ints1,xsorted_infs);

length = max(size(xsorted_infs,1),size(xsorted_infs,2))-2;

inds = zeros(1,length);

k = - nint;
for i=1:length
    k = k+h(i);
    inds(indices(i)) = abs(k);
end

