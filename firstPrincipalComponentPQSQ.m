function [V,U,C,explainedVariance,MDS] = firstPrincipalComponentPQSQ(x, intervals, potential_function_handle, varargin)

verbose=0;
eps=0.001;
initiatilization = 0; % 0 - PC1, -1 - random vector, >0 - ith data vector
meanDefined = 0;

for i=1:length(varargin)
    if strcmp(varargin{i},'verbose')
        verbose = varargin{i+1};
    end
    if strcmp(varargin{i},'eps')
        eps = varargin{i+1};
    end
    if strcmp(varargin{i},'init')
        initiatilization = varargin{i+1};
    end
    if strcmp(varargin{i},'mean')
        meanDefined = 1;
        meanVector = varargin{i+1};
    end
end

intervals_inf = zeros(size(x,2),size(intervals,2)+1);

for k=1:size(x,2) 
    intervals_inf(k,1:size(intervals,2)) = intervals(k,:);
    intervals_inf(k,size(intervals,2)+1) = Inf;
end


[A,B] = computeABcoefficients(intervals, potential_function_handle);

% if meanDefined
% C = meanVector;
% else
% C = PQSQ_Mean(x,intervals,potential_function_handle);
% end
V = zeros(1,size(x,2));
V1 = zeros(1,size(x,2));
U = zeros(size(x,1),1);
U1 = U;
N = size(x,1);
%mn = C;
stdev = std(x);

RS = zeros(size(x,1),size(x,2));



%%%%%%%%%%%%%%%%%%%%%%%%%%
% Center data matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xc = zeros(size(x,1),size(x,2));
% for i=1:size(x,1) 
%     Xc(i,:)= x(i,:)-mn; 
% end;

Xc = x;



XU = 0;
SU = 0;
SU2 = 0;

%%%%%%%%%%%%%%%%%%%%%
% initialize V and U
%%%%%%%%%%%%%%%%%%%%%
%for k=1:size(x,2)
%    V(k) = rand()*stdev(k);
%end
%V = firstPrincipalComponentStandardL2(x);
%V =[0.4451   -0.2209    0.4951    0.4817];

if initiatilization==-1
    for k=1:size(x,2)
        V(k) = rand()*stdev(k);
    end
end

if initiatilization==0
    pc = pca(Xc);
    pc1 = pc(:,1);
    for i=1:size(pc1,1)
        V(i) = pc1(i);
    end
end

if initiatilization>0
    for i=1:size(x,2)
        V(i) = Xc(initiatilization,i);
    end
end

for i=1:size(x,1)
    U(i)=sum(Xc(i,:).*V)/sum(V.*V);
end

count = 1;

V = V/norm(V);

while(count<1000)

U1 = U;

for k=1:size(x,2)
    %inds = splitPointIntervals(Xc(:,k)-V(k)*U(:),intervals(k,:));
    %[x1cinf,indices,intervals_symm] = prepareForSplitPointIntervalsFast(Xc(:,k)-V(k)*U(:),intervals(k,:));
    %inds = splitPointIntervalsFast(x1cinf,indices,0,intervals_symm);
    inds = splitPointIntervalsFast1(Xc(:,k)-V(k)*U(:),intervals_inf(k,:));
    RS(:,k) = inds(:);
end

for i=1:size(x,1)
    XU=0;
    SU2=0;
    for k=1:size(x,2)
        XU=XU+A(k,RS(i,k))*Xc(i,k)*V(k);
        SU2=SU2+A(k,RS(i,k))*V(k)*V(k);
    end
     if SU2==0
         U(i)=0;
     else
         U(i)=XU/SU2;
     end
end


for k=1:size(x,2)
    
XV = 0;
SV2 = 0;
    
    for i=1:size(x,1)
        XV=XV+A(k,RS(i,k))*Xc(i,k)*U(i);
        SV2=SV2+A(k,RS(i,k))*U(i)*U(i);
    end
    if SV2==0
    V1(k) = 0;
    else
    V1(k) = XV/SV2;
    end
end

V1 = V1/norm(V1);

delta = (norm(V-V1))/norm(V1);

if(delta<eps)
    V = V1;
    break;
end

V=V1;

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

if(V(1)<0) 
    V=-V;
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
