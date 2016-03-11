function y = PQSQ(x, intervals, potential_function_handle)


for k=1:size(x,1)
    
R = max(intervals(:,:));

a = 0;
b = 0;

        for i=1:size(intervals,2)-1
        if((intervals(1,i)<=abs(x(k)))&(abs(x(k))<intervals(1,i+1)))
                xk = intervals(1,i);
                xk1 = intervals(1,i+1);
                a = (potential_function_handle(xk)-potential_function_handle(xk1))/(xk*xk-xk1*xk1);
                b = (potential_function_handle(xk1)*xk*xk-potential_function_handle(xk)*xk1*xk1)/(xk*xk-xk1*xk1);
         end
    end
    
    if(abs(x(k))>=intervals(1,size(intervals,2)))
        a = 0;
        b = potential_function_handle(R);
        %display(sprintf('x(%f)>R(%f)!',x(k),R));
    end
    
    y(k) = b+a*x(k)*x(k);
end

end