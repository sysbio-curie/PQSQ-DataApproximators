function c = PQSQ_Mean_Prime(x, intervals, potential_function_handle, mn, varargin)
%PQSQ_Mean calculates PQSQ central point: point which provide the minimum
%of PQSQ error sum for x-c.

    verbose=0;
    eps=0.001;
    
    m = size(x, 2);
    c = zeros(1, m);

    %Form complete intervals by adding the last column of positive infinity
    intervals_inf = [intervals, Inf(m,1)];


    for i=1:length(varargin)
        if strcmp(varargin{i},'verbose')
            verbose = varargin{i+1};
        end
        if strcmp(varargin{i},'eps')
            eps = varargin{i+1};
        end
    end
    
    [A,~] = computeABcoefficients(intervals, potential_function_handle);
    
    for k=1:m
        
        if(verbose)
            display(sprintf('Coordinate %i:',k));
        end
        
        m = mn(k);
        
        count=0;
        
        while(count<100)
            
            m0=m;
            inds = splitPointIntervalsFast1(abs(x(:,k)-m),intervals_inf(k,:));
            as = A(k,inds)';
            x2 = sum(as);
            x1 = sum(as.*x(:,k));
            if x2~=0
                m=x1/x2;
            else
                m=0;
            end
            
            count=count+1;
            
%             if(verbose)
%                 dist_PQSQ=0; dist_fx=0;
%                 for i=1:size(x,1)
%                     pqsq = PQSQ(x(i,k)-m,intervals(k,:),potential_function_handle);
%                     fd = potential_function_handle(x(i,k)-m);
%                     %display(sprintf('%f %f',pqsq,fd));
%                     dist_PQSQ=dist_PQSQ+abs(pqsq);
%                     dist_fx=dist_fx+abs(fd);
%                 end
%                 dist_PQSQ=dist_PQSQ/size(x,1);
%                 dist_fx=dist_fx/size(x,1);
%                 
%                 display(sprintf('%i: m=%f Error in PQSQ=%f, Error in f(x)=%f',count,m,dist_PQSQ,dist_fx));
%             end
            
            delta = abs(m-m0)/abs(m+0.001);
            %display(sprintf('%i: Delta=%f',count,delta));
            if(delta<eps)
                c(k) = m;
                break;
            end;
            
            %display(sprintf('%i: ',count));
            
        end;
        
%         if(verbose)
%             dm=0;
%             dmean = 0;
%             md = median(x(:,k));
%             meand = mean(x(:,k));
%             for i=1:size(x,1)
%                 dm=dm+sum(abs(md-x(i,k)));
%                 dmean=dmean+sum(abs(meand-x(i,k)));
%             end
%             dm=dm/size(x,1);
%             dmean=dmean/size(x,1);
%             display(sprintf('Reference: L1 distance to median=%f, L1 distance to mean=%f',dm,dmean));
%         end
        
        %y(:,k)=m(:);
        
    end;
    
end

