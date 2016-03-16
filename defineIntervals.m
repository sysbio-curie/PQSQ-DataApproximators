function intervals = defineIntervals(x, number_of_intervals)
%defineIntervals define "uniform" intervals for data matrix x and specified
%number of intervals.
% x is n-by-m matrix.
% number_of_intervals is required number of intervals.
%Output argument intervals is matrix m-by-number_of_intervals. Each row
%   contains number_of_intervals values of thresholds for intervals.

    p=number_of_intervals-1;

    %Calculate characteristic distance for all coordinates
    characteristic_distance = (max(x)-min(x))';

    %Calculate j^2/p^2
    row = ((0:p)/p).^2;
    
    %Now intervals is the product of row and characteristic_distance:
    intervals = characteristic_distance * row;
    
    
%     for i=1:size(x,2)
% 
%         %characteristic_distance = mad(x(:,i))*6;
%         characteristic_distance = max(x(:,i))-min(x(:,i));
%         %characteristic_distance = characteristic_distance*0.1;
% 
%         delta = sqrt(characteristic_distance)/number_of_intervals;
% 
%         s = 0:delta:sqrt(characteristic_distance);
% 
%         intervals(i,:) = s(:);
% 
%     end;
% 
%     intervals(:) = intervals(:).*intervals(:);

    %if intervals(i,size(intervals,2))<max(x(:,i))*2
    %intervals(i,size(intervals,2)) = max(x(:,i))*2;
    %end


end