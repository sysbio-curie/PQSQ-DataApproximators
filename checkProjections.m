k = 1; 

clear distL1; clear distL2; clear distPQSQ; 
clear norm1; clear norm2;
d = U(k)-100:0.01:U(k)+100;
distL1 = zeros(size(d,2),1); distL2 = zeros(size(d,2),1); distPQSQ = zeros(size(d,2),1);
norm1 = zeros(size(d,2),1); norm2 = zeros(size(d,2),1);
for i=1:size(d,2) [distPQSQ(i),norm1(i),norm2(i)]=PQSQ_Norm(x(k,:)-C-d(i)*V',intervals,@L1); distL2(i)=norm(x(k,:)-C-d(i)*V'); distL1(i)=sum(abs(x(k,:)-C-d(i)*V')); end;
plot(d,distPQSQ,'b-','LineWidth',3); hold on; 

plot(d,norm1,'b-','LineWidth',1); hold on; 
plot(d,norm2,'b-','LineWidth',1); hold on; 

plot(d,distL2,'g-','LineWidth',2);  
plot(d,distL1,'m-','LineWidth',2);  
plot([U(k) U(k)],[0 max(distPQSQ)],'r--','LineWidth',3); 

 dPQSQ = d(find(distPQSQ == min(distPQSQ))) 
 dL2 = d(find(distL2 == min(distL2)))
 dL1 = d(find(distL1 == min(distL1)))
 U(k)
 plot([dPQSQ dPQSQ],[0 -1],'b--','LineWidth',2); 
 plot([dL2 dL2],[0 -1],'g--','LineWidth',2); 
 plot([dL1 dL1],[0 -1],'m--','LineWidth',2);
 plot([min(d) max(d)],[0 0 ],'k-');