function [V,U] = onePC_PQSQR(x, potentialFunction, varargin)
%onePC_PQSQR calculates the first principal component for specified data
%matrix x and PQSQ potential function
%   x is n-by-m matrix with n observations and m attributes.
%   potentialFunction is structure which define PQSQ potential function
%       (see definePotentialFunction.m for details)
% Name, Value optionas:
% 'eps', positive number is tolerance level for component calculation.
% 'init', -1/0/positive number is method of component initialization:
%               -1 means uniform random generation in interval 0, stdev for
%                   each coordiate
%               0 means initialization by the first L2 principal component
%                   for data matrix x.
%               positive number is the number of data point to use it as
%                   initialization.
%Output arguments
%V is m-by-1 vector of found principal component.
%U is n-by-1 vector of length of projections onto found principal
%   component.

    eps=0.001;
    initiatilization = 0; % 0 - PC1, -1 - random vector, >0 - ith data vector

    for i=1:length(varargin)
        if strcmp(varargin{i},'eps')
            eps = varargin{i+1};
        elseif strcmp(varargin{i},'init')
            initiatilization = varargin{i+1};
        end
    end

    %define dimensions of vectors and matrices
    m = size(x, 2);

    %%%%%%%%%%%%%%%%%%%%%
    % initialize V and U
    %%%%%%%%%%%%%%%%%%%%%
    if initiatilization==-1
        stdev = std(x)';
        V = rand(m,1).*stdev;
    elseif initiatilization==0
        [V,~] = eigs(x'*x,1);
    elseif initiatilization>0
        V=x(initiatilization,:)';
    end

    %Renormalize V
    V = V/norm(V);
    %Initialize U
    U = x*V;
    count = 1;
    while(count<1000)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %U recalculation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate the indices of intervals
        inds = identifyIntervals(sum((x-U*V').^2,2),potentialFunction.sqint);
        
        %Identify nonzero coefficients a
        ax = potentialFunction.A(inds);
        inds = ax == 0;
        
        %Calculate intermediate results
        XV = (x*V)/(V'*V);
        
        %Finilize
        U = XV;
        U(inds) = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Recalculate V. Remember old value
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        V1=V;

        %calculate indices of intervals and new principal component
        %calculate the indices of intervals
        inds = identifyIntervals(sum((x-U*V').^2,2),potentialFunction.sqint);
        ax = potentialFunction.A(inds)';
        SAU = sum(ax.*(U.^2));
        SXU = sum(bsxfun(@times,x,ax.*U));
        inds = SXU == 0;
        V(inds) = 0;
        V(~inds) = SXU(~inds)./SAU(~inds);
        
        %renormalize V
        V = V/norm(V);
        
        delta = (norm(V-V1));
        
        if(delta<eps)
            break;
        end
        
        count = count+1;
    end
    
    %Standartize direction
    if(V(1)<0)
        V=-V;
        U=-U;
    end
    
    %Final projection recalculation
    inds = identifyIntervals(sum((x-U*V').^2,2),potentialFunction.sqint);

    %Identify nonzero coefficients a
    ax = potentialFunction.A(inds)';
    inds = ax == 0;

    %Calculate intermediate results
    XV = (x*V)/(V'*V);

    %Finilize
    U = XV;
    U(inds) = 0;
end