function intervals = defineIntervals(x, number_of_intervals)

    for i=1:size(x,2)

   %characteristic_distance = mad(x(:,i))*6;
   characteristic_distance = max(x(:,i))-min(x(:,i));
   %characteristic_distance = characteristic_distance*0.1;
    
    delta = sqrt(characteristic_distance)/number_of_intervals;
    
    s = 0:delta:sqrt(characteristic_distance);
    
    intervals(i,:) = s(:);
    
    end;
    
    intervals(:) = intervals(:).*intervals(:);

    %if intervals(i,size(intervals,2))<max(x(:,i))*2
    %intervals(i,size(intervals,2)) = max(x(:,i))*2;
    %end
    
    
end