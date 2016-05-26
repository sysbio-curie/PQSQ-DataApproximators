function intervals = PQSQR_defineIntervals(x, number_of_intervals, delta)
%PQSQR_defineIntervals define "uniform in squares" intervals for data
%matrix x and specified number of intervals.
%   x is n-by-m matrix.
%   number_of_intervals is required number of intervals.
%   delta is coefficient of shrinkage which is greater than 0 ang not
%       greater than 1 (1 means without trimming).
%Output argument intervals is row vector 1-by-number_of_intervals which
%   contains number_of_intervals values of thresholds for intervals.

    p=number_of_intervals-1;
    
    if nargin<3
        delta = 1;
    end

    %Calculate distances between all points
    d=sum(x.^2,2);
    D = bsxfun(@plus,d,d')-2*(x*x');
    
    %Maximal distance is characteristic distance
    characteristic_distance = sqrt(max(D(:)))*delta;
    
    %Calculate j^2/p^2
    row = ((0:p)/p).^2;
    
    %Now intervals is the product of row and characteristic_distance:
    intervals = characteristic_distance * row;
end