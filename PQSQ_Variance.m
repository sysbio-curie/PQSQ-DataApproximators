function var = PQSQ_Variance(x, intervals, potential_function_handle)

var = zeros(1,size(x,2));

PQSQmean = PQSQ_Mean(x,  intervals, potential_function_handle);


for k=1:size(x,2)

xx = PQSQ(x(:,k)-PQSQmean(k),intervals(k,:),potential_function_handle)/size(x,1);

for i=1:size(x,1)
var(k) = var(k)+xx(i);
end

var(k)/size(x,1);

end

end