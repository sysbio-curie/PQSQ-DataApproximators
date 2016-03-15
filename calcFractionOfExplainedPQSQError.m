function PQSQerr = calcFractionOfExplainedPQSQError(x,V,U,C,intervals,potential_function_handle)

C1 = repmat(C,size(x,1),1);
xproj = C1+U*V;
errunexplained = calcAveragePQSQError(x,xproj,intervals,potential_function_handle);
errtotal =calcAveragePQSQError(x,C1,intervals,potential_function_handle); 

PQSQerr = (errtotal-errunexplained)/errtotal;

end

