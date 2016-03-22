function U = optimizePQSQ1D(x,V,Uguess,C,intervals,potential_function_handle)
    
    m = size(x,2);
    Xc = x - C;
    [A,~] = computeABcoefficients(intervals, potential_function_handle);
    intervals_inf = [intervals, Inf(m,1)];

for i=1:m
    uk = Xc(i)/V(i);
    
    for j=1:size(x,2)
        inds(j) = splitPointIntervalsFast1(Xc(j)-V(j)*uk,intervals_inf(j,:));
    end
    XU = 0;
    SU = 0;
    for j=1:m
            XU = XU + A(j,inds(j))*(Xc(j))*V(j);
            SU = SU + A(j,inds(j))*V(j)*V(j);
    end;
    ukopt(i) = XU/SU;
end

for i=1:m
    err(i) = PQSQ_Norm(Xc-ukopt(i)*V',intervals,@L1);
end

U1 = ukopt(err==min(err));
U = U1(1);

    
end

