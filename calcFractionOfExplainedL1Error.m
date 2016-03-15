function L1av = calcFractionOfExplainedL1Error(x,V,U,C)
%Returns explained average absolute error to a basis

C1 = repmat(C,size(x,1),1);
xproj = C1+U*V;
errunexplained = calcAverageL1Distance(x,xproj);
errtotal =calcAverageL1Distance(x,C1); 

L1av = (errtotal-errunexplained)/errtotal;

end

