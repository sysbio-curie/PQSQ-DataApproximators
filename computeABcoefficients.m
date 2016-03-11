function [A,B] = computeABcoefficients(intervals, potential_function_handle)

A = zeros(size(intervals,1),size(intervals,2));
B = zeros(size(intervals,1),size(intervals,2));

for k=1:size(intervals,1)

    for i=1:size(intervals,2)-1
                xk = intervals(k,i);
                xk1 = intervals(k,i+1);
                A(k,i) = (potential_function_handle(xk)-potential_function_handle(xk1))/(xk*xk-xk1*xk1);
                B(k,i) = (potential_function_handle(xk1)*xk*xk-potential_function_handle(xk)*xk1*xk1)/(xk*xk-xk1*xk1);
    end
    
    A(k,size(intervals,2)) = 0;
    B(k,size(intervals,2)) = potential_function_handle(intervals(k,size(intervals,2)));
    
end    
    
end