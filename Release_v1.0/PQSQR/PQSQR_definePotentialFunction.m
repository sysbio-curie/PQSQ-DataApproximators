function potentialFunction = PQSQR_definePotentialFunction( x, number_of_intervals, potential_function_handle, delta )
%PQSQR_definePotentialFunction defines "uniform in square" intervals for
%data matrix x and specified number_of_intervals.
%   x is n-by-m matrix.
%   number_of_intervals is required number of intervals OR row array with
%       intervals boundaries which are specified by user.
%   potential_function_handle is function handler for coefficients
%       calculation.
%   delta is coefficient of shrinkage which is greater than 0 and not
%       greater than 1 (1 means without trimming).
%Output argument potentialFunction is structure with three fields:
%   intervals is row vector 1-by-number_of_intervals which contains
%       number_of_intervals values of thresholds for intervals and one
%       additional value Inf
%   sqint is row vector 1-by-number_of_intervals with squared intervals
%       boundaries.
%   A and B are the row vectors 1-by-number_of_intervals with quadratic
%       functions coefficients

    if nargin<4 
        delta = 1;
    end
    
    if isscalar(number_of_intervals)
        intervals = PQSQR_defineIntervals(x, number_of_intervals, delta);
    else
        intervals = (number_of_intervals(:))';
    end
    
    potentialFunction.intervals = [intervals, Inf(1)];
    potentialFunction.sqint = potentialFunction.intervals.^2;
    [potentialFunction.A,potentialFunction.B] = ...
        PQSQR_computeABcoefficients(intervals, potential_function_handle);
end