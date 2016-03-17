function [U,RS] = optimizeProjections(Xcentered,V,U,A,intervals_inf,potential_function_handle,MAXNUMBER_OF_ITERATIONS,EPS)

[n, m] = size(Xcentered);
RS = zeros(n, m);
XU=U;
SU=U;
    
 ktest = 15;
 d = min(U)-1.2:0.05:max(U)+1.2;
 for i=1:size(d,2) 
     distPQSQ(i)=PQSQ_Norm(Xcentered(ktest,:)-d(i)*V',intervals_inf,@L1); 
 end;
 trueUtest = d(find(distPQSQ == min(distPQSQ)));
 Utest = U(ktest);

% d = min(U):0.1:max(U);
% for ktest = 1:size(Xcentered,1)
%     %display(sprintf('%i',ktest));
%  for i=1:size(d,2) 
%      distPQSQ(i)=PQSQ_Norm(Xcentered(ktest,:)-d(i)*V',intervals_inf,@L1); 
%  end;
%  trueUtest = d(find(distPQSQ == min(distPQSQ)));
%  U(ktest) = trueUtest;
%  Utest = U(ktest);
% end
% 
%         for k=1:m
%             inds = splitPointIntervalsFast1(Xcentered(:,k)-V(k)*U(:),intervals_inf(k,:));
%             RS(:,k) = inds(:);
%         end
%count = MAXNUMBER_OF_ITERATIONS+1;

count = 1;

while(count<=MAXNUMBER_OF_ITERATIONS)

        U1 = U;

        U = U+rand(n,1)*0.1;
        
        for k=1:m
            inds = splitPointIntervalsFast1(Xcentered(:,k)-V(k)*U(:),intervals_inf(k,:));
            RS(:,k) = inds(:);
        end

        
        XU(:)=0;
        SU(:)=0;
        for k=1:m
            XU = XU + A(k,RS(:,k))'.*Xcentered(:,k)*V(k);
            SU = SU + A(k,RS(:,k))'*V(k)*V(k);
        end
        ind = SU==0;
        U(ind)=0;
        U(~ind)=XU(~ind)./SU(~ind);
        
        delta = (norm(U-U1))/norm(U);
        
        if(delta<EPS)
            break;
        end

count=count+1;        
        
end;

display(sprintf('trueUtest=%f,Utest=%f,U(ktest)=%f',trueUtest,Utest,U(ktest)));

end

