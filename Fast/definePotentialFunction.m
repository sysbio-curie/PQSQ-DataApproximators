function potentialFunction = definePotentialFunction( x, number_of_intervals, potential_function_handle, delta )
%definePotentialFunction defines "uniform" intervals for data matrix x and
%specified number_of_intervals.
%   x is n-by-m matrix.
%   number_of_intervals is required number of intervals.
%   potential_function_handle is function handler for coefficients
%       calculation.
%   delta is coefficient of shrinkage which is greater than 0 ang not
%       greater than 1 (1 means without trimming).
%Output argument potentialFunction is structure with three fields:
%   intervals is matrix m-by-number_of_intervals. Each row contains
%       number_of_intervals values of thresholds for intervals and one
%       additional value Inf
%   A and B are the m-by-number_of_intervals matrices with quadratic
%       functions coefficients

    p=number_of_intervals-1;
    
    if nargin<4 
        delta = 1;
    end

    %Calculate characteristic distance for all coordinates
    characteristic_distance = (max(x)-min(x))'*delta;

    %Calculate j^2/p^2
    row = ((0:p)/p).^2;
    
    %Now intervals is the product of row and characteristic_distance:
    intervals = characteristic_distance * row;
    
    potentialFunction.intervals = [intervals, Inf(size(x,2),1)];
    [potentialFunction.A,potentialFunction.B] = computeABcoefficients_fast(intervals, potential_function_handle);
    
    


end