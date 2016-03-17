function [norm,norm1,norm2]=PQSQ_Norm(vector, intervals, potential_function_handle)

%[A,B] = computeABcoefficients(intervals, potential_function_handle);

norm = 0;

for k=1:max(size(vector,1),size(vector,2))
    normk = PQSQ(vector(k),intervals(k,:),potential_function_handle);
    if k==1 norm1 = normk; end
    if k==2 norm2 = normk; end
    norm = norm + normk;
end