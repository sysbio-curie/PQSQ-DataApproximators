function inds = splitPointIntervals(x,intervals)

dim = max(size(x,1),size(x,2));
inds = zeros(1,dim);

for s=1:dim

k=0;    

val = abs(x(s));


   for i=1:size(intervals,2)-1
        if((intervals(i)<=val)&(val<intervals(i+1)))
            k=i;
            inds(s) = k;
            break; 
        end
   end
   
    if(val>=intervals(size(intervals,2)))
       k = size(intervals,2);
       inds(s) = k;
    end
    
    %if k==0
    %    display('');
    %end
    
end

end
