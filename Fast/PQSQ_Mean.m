function c = PQSQ_Mean(x, potentialFunction, eps)
%PQSQ_Mean calculates PQSQ central point: point which provide the minimum
%of PQSQ error sum for x-c.
%   x is n-by-m matrix of data
%   potentialFunction is structure which define PQSQ potential function
%       (see definePotentialFunction.m for details)
%   eps is tolerance level for iterations convergence

    if nargin<3
        eps = 0.001;
    end

    m = size(x, 2);
    
    %Calculate results for three 'typical' initial points
    center_zero = PQSQ_Mean_Prime(x,potentialFunction,zeros(1,m),eps);
    center_mean = PQSQ_Mean_Prime(x,potentialFunction,mean(x),eps);
    center_median = PQSQ_Mean_Prime(x,potentialFunction,median(x),eps);

    %Calculate norms for comparison
    norm = [PQSQ_Norm(bsxfun(@minus,x,center_zero), potentialFunction,1);...
        PQSQ_Norm(bsxfun(@minus,x,center_mean), potentialFunction,1);...
        PQSQ_Norm(bsxfun(@minus,x,center_median), potentialFunction,1)];
    
    %Select the minimal values
    [~,ind] = min(norm);
    
    %Form output vector
    c = center_zero;
    ii = ind==2;
    c(ii) = center_mean(ii);
    ii = ind==3;
    c(ii) = center_median(ii);
end

