function [V,U,C,explainedVariance,MDS] = firstPrincipalComponentPQSQ(x, intervals, potential_function_handle, varargin)
%firstPrincipalComponentPQSQ calculates the first principal component for
%specified data matrix x, intervals and potential_function_handle.
% x is n-by-m matrix with n observations and m attributes.
% intervals is m-by-(p+1) matrix with (p+1) boundaries of p intervals for
%   each dimension of data in row.
% potential_function_handle is handler for energy potential function.
% Name, Value optionas:
% 'verbose', 1 is indicator of detailed output (debuging).
% 'eps', positive number is tolerance level for component calculation.
% 'init', -1/0/positive number is method of component initialization:
%               -1 means uniform random generation in interval 0, stdev for
%                   each coordiate
%               0 means initialization by the first L2 principal component
%                   for data matrix x.
%               positive number is the number of data point to use it as
%                   initialization.
%
%Output arguments
%V is m-by-1 vector of found principal component.
%U is n-by-1 vector of length of projections onto found principal
%   component.

    verbose=0;
    eps=0.001;
    initiatilization = 0; % 0 - PC1, -1 - random vector, >0 - ith data vector
%     meanDefined = 0;

    for i=1:length(varargin)
        if strcmp(varargin{i},'verbose')
            verbose = varargin{i+1};
        elseif strcmp(varargin{i},'eps')
            eps = varargin{i+1};
        elseif strcmp(varargin{i},'init')
            initiatilization = varargin{i+1};
        end
%         if strcmp(varargin{i},'mean')
%             meanDefined = 1;
%             meanVector = varargin{i+1};
%         end
    end

    %define dimensions of vectors and matrices
    [n, m] = size(x);
    
    %Form complete intervals by adding the last column of positive infinity
    intervals_inf = [intervals, Inf(m,1)];

% intervals_inf = zeros(size(x,2),size(intervals,2)+1);
% 
% for k=1:size(x,2) 
%     intervals_inf(k,1:size(intervals,2)) = intervals(k,:);
%     intervals_inf(k,size(intervals,2)+1) = Inf;
% end


    [A,~] = computeABcoefficients(intervals, potential_function_handle);
%[A,B] = computeABcoefficients(intervals, potential_function_handle);

% if meanDefined
% C = meanVector;
% else
% C = PQSQ_Mean(x,intervals,potential_function_handle);
% end
    %Create auxiliary arrays
%V is m-by-1 vector of found principal component.
%U is 1-by-n vector of length of projections onto found principal
    RS = zeros(n, m);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Center data matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xc = zeros(size(x,1),size(x,2));
% for i=1:size(x,1) 
%     Xc(i,:)= x(i,:)-mn; 
% end;

    Xc = x;

    %%%%%%%%%%%%%%%%%%%%%
    % initialize V and U
    %%%%%%%%%%%%%%%%%%%%%
    if initiatilization==-1
        stdev = std(x)';
        V = rand(m,1).*stdev;
    %     for k=1:size(x,2)
    %         V(k) = rand()*stdev(k);
    %     end
    elseif initiatilization==0
        [V,~] = eigs(Xc'*Xc,1);
    elseif initiatilization>0
        V=Xc(initiatilization,:)';
    end

    %Renormalize V
    V = V/norm(V);
    %Initialize U
    U = Xc*V;
       
    
    count = 1;
    
    XU=U;
    SU=U;
    XV = V;
    SV = V;
    
    while(count<1000)
        
        %Copy old vector U
        U1 = U;
        
        %calculate the matrix of indeces of intervals
        for k=1:m
            %inds = splitPointIntervals(Xc(:,k)-V(k)*U(:),intervals(k,:));
            %[x1cinf,indices,intervals_symm] = prepareForSplitPointIntervalsFast(Xc(:,k)-V(k)*U(:),intervals(k,:));
            %inds = splitPointIntervalsFast(x1cinf,indices,0,intervals_symm);
            inds = splitPointIntervalsFast1(Xc(:,k)-V(k)*U(:),intervals_inf(k,:));
            RS(:,k) = inds(:);
        end
        
% This fragment approximately 40 times slower than following new fragment        
%         for i=1:size(x,1)
%             XU=0;
%             SU2=0;
%             for k=1:size(x,2)
%                 XU=XU+A(k,RS(i,k))*Xc(i,k)*V(k);
%                 SU2=SU2+A(k,RS(i,k))*V(k)*V(k);
%             end
%             if SU2==0
%                 U(i)=0;
%             else
%                 U(i)=XU/SU2;
%             end
%         end

        %Initialize arrays XU and SU
        
        XU(:)=0;
        SU(:)=0;
        for k=1:m
            XU = XU + A(k,RS(:,k))'.*Xc(:,k)*V(k);
            SU = SU + A(k,RS(:,k))'*V(k)*V(k);
        end
        ind = SU==0;
        U(ind)=0;
        U(~ind)=XU(~ind)./SU(~ind);
        
        V1=V;
% This fragment is approximately 25-30 times slower than following new
%         V1=V;
%         for k=1:size(x,2)
%             XV = 0;
%             SV2 = 0;
%             for i=1:size(x,1)
%                 XV=XV+A(k,RS(i,k))*Xc(i,k)*U(i);
%                 SV2=SV2+A(k,RS(i,k))*U(i)*U(i);
%             end
%             if SV2==0
%                 V1(k) = 0;
%             else
%                 V1(k) = XV/SV2;
%             end
%         end

        for k=1:size(x,2)
            XV(k) = sum(A(k,RS(:,k))'.*Xc(:,k).*U(:));
            SV(k) = sum(A(k,RS(:,k))'.*U(:).*U(:));
        end
        ind = SV==0;
        V1(ind)=0;
        V1(~ind)=XV(~ind)./SV(~ind);
        
        
        V1 = V1/norm(V1);
        
        delta = (norm(V-V1));
        
        V=V1;
        if(delta<eps)
            break;
        end
        
        
        count = count+1;
        if verbose
            PQSQ_mds = 0;
            deltaU = (norm(U-U1))/norm(U1);
            for i=1:size(x,1)
                diff = PQSQ_Norm(Xc(i,:)-U(i)*V,intervals,potential_function_handle);
                PQSQ_mds=PQSQ_mds+diff;
            end
            PQSQ_mds = PQSQ_mds/N;
            MDS(count) = PQSQ_mds;
            
            display(sprintf('%i: DeltaV=%f, DeltaU=%f, PQSQ_MDS=%f', count, delta, deltaU, PQSQ_mds));
            
        end
    end
    
    %Standartize direction
    if(V(1)<0)
        V=-V;
        U=-U;
    end
    
    if verbose
        PQSQ_mds = 0;
        for i=1:size(x,1)
            diff = PQSQ_Norm(Xc(i,:)-U(i)*V,intervals,potential_function_handle);
            PQSQ_mds=PQSQ_mds+diff;
        end
        PQSQ_mds = PQSQ_mds/N;
        totalVariance = sum(PQSQ_Variance(x,intervals,potential_function_handle));
        explainedUariance = totalVariance - PQSQ_mds;
        
        display(sprintf('Fraction of explained PQSQ variance: %f', explainedUariance/totalVariance));
        
    end
