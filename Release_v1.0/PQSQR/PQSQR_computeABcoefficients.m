function [A,B] = PQSQR_computeABcoefficients(intervals, potential_function_handle)
%PQSQR_computeABcoefficients calculates the coefficients a and b for
%quadratic fragments of potential function.
%   intervals is the 1-by-K matrix of intervals' boudaries without final
%       infinit boundary.
%   potential_function_handle is a handle of majorant function.

    %Get dimensions of intervals
    p = size(intervals,2);

    %Preallocate memory
    A = zeros(1,p);
    B = zeros(1,p);

    %Calculate value of function all boundaries
    pxk = potential_function_handle(intervals);
    sxk = intervals.^2;

    A(1:p-1) = (pxk(1:p-1)-pxk(2:p))./(sxk(1:p-1)-sxk(2:p));
    B(1:p-1) = (pxk(2:p).*sxk(1:p-1)-pxk(1:p-1).*sxk(2:p))./...
        (sxk(1:p-1)-sxk(2:p));
    B(p) = pxk(p);
end