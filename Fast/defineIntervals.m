function intervals = defineIntervals(x, number_of_intervals, delta)
%defineIntervals define "uniform" intervals for data matrix x and specified
%number of intervals.
%   x is n-by-m matrix.
%   number_of_intervals is required number of intervals.
%   delta is coefficient of shrinkage which is greater than 0 ang not
%       greater than 1 (1 means without trimming).
%Output argument intervals is matrix m-by-number_of_intervals. Each row
%   contains number_of_intervals values of thresholds for intervals.

    p=number_of_intervals-1;
    
    if nargin<3
        delta = 1;
    end

    %Calculate characteristic distance for all coordinates
    characteristic_distance = (max(x)-min(x))'*delta;

    %Calculate j^2/p^2
    row = ((0:p)/p).^2;
    
    %Now intervals is the product of row and characteristic_distance:
    intervals = characteristic_distance * row;
end