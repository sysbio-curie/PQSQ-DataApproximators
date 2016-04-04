function c = PQSQ_Mean(x, intervals, potential_function_handle, varargin)
%PQSQ_Mean calculates PQSQ central point: point which provide the minimum
%of PQSQ error sum for x-c.

    m = size(x, 2);
    
    %Calculate results for three 'typical' initial points
    center_zero = PQSQ_Mean_Prime(x,intervals,potential_function_handle,zeros(1,m),varargin);
    center_mean = PQSQ_Mean_Prime(x,intervals,potential_function_handle,mean(x),varargin);
    center_median = PQSQ_Mean_Prime(x,intervals,potential_function_handle,median(x),varargin);

    %Calculate norms for comparison
    norm = [sum(PQSQ_Norm(bsxfun(@minus,x,center_zero), intervals, potential_function_handle));...
        sum(PQSQ_Norm(bsxfun(@minus,x,center_mean), intervals, potential_function_handle));...
        sum(PQSQ_Norm(bsxfun(@minus,x,center_median), intervals, potential_function_handle))];
    
    %Select the minimal values
    [~,ind] = min(norm);
    
    %Form output vector
    c = center_zero;
    ii = ind==1;
    c(ii) = center_mean(ii);
    ii = ind==2;
    c(ii) = center_median(ii);
end

