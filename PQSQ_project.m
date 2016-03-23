function b = PQSQ_project(x, v, intervals, funct, initVector)
%PQSQ_project search vector b which provide the minimal
%sum(PQSQ_Norm(x-bv')).
%   x is n-by-m matrix of data points to project
%   v is m-by-1 vector of line direction
%   intervals is matrix of intervals for PQSQ
%   funct is handler of majorant function.
%
%   b is n-by-1 vector of length of projections

    %Calculate initial value of b
    b = initVector;
    XU = b;
    SU = b;
    
    %get sizes
    [n, m] = size(x);

    %Form complete intervals by adding the last column of positive infinity
    intervals_inf = [intervals, Inf(m,1)];

    %Create auxiliary arrays
    [A,~] = computeABcoefficients(intervals, funct);

    %Loop till convergence or at most 1000 iterations
    for qq=1:1000
        %save old version of b
        b1 = b;

        %Initialize arrays XU and SU
        XU(:)=0;
        SU(:)=0;

        %calculate the matrix of indeces of intervals
        for k=1:m
            inds = splitPointIntervalsFast1(x(:,k)-v(k)*b,intervals_inf(k,:));
            XU = XU + A(k,inds)'.*x(:,k)*v(k);
            SU = SU + A(k,inds)'*v(k)*v(k);
        end
        ind = SU==0;
        b(ind)=0;
        b(~ind)=XU(~ind)./SU(~ind);
        
        %check the convergence
        acc = max(abs(b-b1));
        if acc<0.001
            break;
        end
    end

end