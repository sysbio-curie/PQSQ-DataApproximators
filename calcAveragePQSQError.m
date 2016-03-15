function PQSQavError = calcAveragePQSQError(x,xprojected,intervals,potential_function_handle)
% Calculates average approximation error in PQSQ measure

PQSQavError = 0;
xd = x - xprojected;

for i=1:size(x,1)
    PQSQavError=PQSQavError+PQSQ_Norm(xd(i,:),intervals,potential_function_handle);
end
PQSQavError=PQSQavError/size(x,1);


end

