function ints = makeSymmetricIntervals(intervals)
    n = size(intervals,2);
    for i=1:n
        ints(i) = -intervals(n-i+1);
    end 
    for i=1:n
        ints(n+i) = intervals(i);
    end
end