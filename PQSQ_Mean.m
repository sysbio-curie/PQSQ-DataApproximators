function y = PQSQ_Mean(x, intervals, potential_function_handle,varargin)

verbose=0;
eps=0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare matrix for fast interval indexing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indices = zeros(size(x,1),size(x,2));
% for k=1:size(x,2) 
%     [x1cinf,indices_k,intervals_symm_k] = prepareForSplitPointIntervalsFast(x(:,k),intervals(k,:));
%     Xc_sorted_inf(:,k) = x1cinf;
%     indices(:,k) = indices_k;
%     intervals_symm(k,:) = intervals_symm_k;
% end;

intervals_inf = zeros(size(x,2),size(intervals,2)+1);

for k=1:size(x,2) 
    intervals_inf(k,1:size(intervals,2)) = intervals(k,:);
    intervals_inf(k,size(intervals,2)+1) = Inf;
end


for i=1:length(varargin)
    if strcmp(varargin{i},'verbose')
        verbose = varargin{i+1};
    end
    if strcmp(varargin{i},'eps')
        eps = varargin{i+1};
    end
end
    

    mn = mean(x);
    dim = size(mn,2);
    nint = size(intervals,2);

    
    [A,B] = computeABcoefficients(intervals, potential_function_handle);
    
    for k=1:size(x,2)
    
    if(verbose)
    display(sprintf('Coordinate %i:',k));
    end
        
    m = mn(k);
    %m = min(x(:,k))+rand()*(max(x(:,k))-min(x(:,k)));
    
    count=0;
    
    while(count<100)
        
    m0=m;

    distances = dist(x(:,k),m);
    
    %inds = splitPointIntervals(distances,intervals(k,:));
    %inds = splitPointIntervalsFast(Xc_sorted_inf(:,k),indices(:,k),m,intervals_symm(k,:));
    inds = splitPointIntervalsFast1(distances,intervals_inf(k,:));
    %inds1-inds'
    
    x1=0;
    x2=0;
    
    %for s=1:nint-1
    %    as = A(k,s);
    %    x1 = x1+as*sum(x(inds==s,k));
    %    x2 = x2+as*sum(inds==s);
    %end
    
    for i=1:size(x,1)
        as = A(k,inds(i));
        x1 = x1+as*x(i,k);
        x2 = x2+as;
    end
    
    if x2~=0
    m=x1/x2;
    else
    m=0;
    end

    count=count+1;
    
    if(verbose)
    dist_PQSQ=0; dist_fx=0;
    for i=1:size(x,1)
        pqsq = PQSQ(x(i,k)-m,intervals(k,:),potential_function_handle);
        fd = potential_function_handle(x(i,k)-m);
        %display(sprintf('%f %f',pqsq,fd));
        dist_PQSQ=dist_PQSQ+abs(pqsq);
        dist_fx=dist_fx+abs(fd);
    end
    dist_PQSQ=dist_PQSQ/size(x,1);
    dist_fx=dist_fx/size(x,1);
    
    display(sprintf('%i: m=%f Error in PQSQ=%f, Error in f(x)=%f',count,m,dist_PQSQ,dist_fx));
    end
    
    delta = abs(m-m0)/abs(m+0.001);
    %display(sprintf('%i: Delta=%f',count,delta));
    if(delta<eps)
    y(k) = m;
    break;
    end;
    
    %display(sprintf('%i: ',count));
    
    end;
    
    if(verbose)
    dm=0;
    dmean = 0;
    md = median(x(:,k));
    meand = mean(x(:,k));
    for i=1:size(x,1)
        dm=dm+sum(abs(md-x(i,k)));
        dmean=dmean+sum(abs(meand-x(i,k)));
    end
    dm=dm/size(x,1);
    dmean=dmean/size(x,1);
    display(sprintf('Reference: L1 distance to median=%f, L1 distance to mean=%f',dm,dmean));
    end
    
    %y(:,k)=m(:);
    
    end;
    
end

