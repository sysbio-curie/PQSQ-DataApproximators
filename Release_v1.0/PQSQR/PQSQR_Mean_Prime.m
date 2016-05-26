function c = PQSQR_Mean_Prime(x, potentialFunction, mn, eps)
%PQSQR_Mean_Prime calculates PQSQ central point: point which provide the
%minimum of PQSQ error sum for x-c.
%   x is n-by-m matrix
%   potentialFunction is structure which define PQSQ potential function
%       (see PQSQR_definePotentialFunction.m for details)
%   mn is the initial state for centre.
%   eps is tolerance level for iterations convergence

    c = (mn(:))';
    xx = sum(x.^2,2);
    count=0;
    while(count<100)
        %Calculate distances from data points to current centre
        mn = c;
        cc=sum(c.^2);
        d=xx+cc-2*(x*c');
        %Identify intervals for each points
        inds = identifyIntervals(d,potentialFunction.sqint);
        as = potentialFunction.A(inds)';
        x2 = sum(as);
        x1 = sum(bsxfun(@times,x,as));
        if x2~=0
            c = x1./x2;
        else
            c = 0;
        end
        
        count = count+1;
        delta = norm(c-mn)/(norm(c)+0.001);
        if(delta<eps)
            break;
        end;
    end;
end