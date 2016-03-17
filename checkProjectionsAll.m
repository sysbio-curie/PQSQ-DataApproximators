clear distL1; clear distL2; clear distPQSQ; clear uL2;
d = min(U)-1.2:0.05:max(U)+1.2;
distL1 = zeros(size(d,2),1); 
distL2 = zeros(size(d,2),1); 
distPQSQ = zeros(size(d,2),1);
uL2 = zeros(size(d,2),1);


for k = 1:size(x,1)

display(sprintf('%i',k));
uL2(k) = sum((x(k,:)-C).*V')/sum(V.*V);
    
for i=1:size(d,2) 
    distPQSQ(i)=PQSQ_Norm(x(k,:)-C-d(i)*V',intervals,@L1); 
    distL2(i)=norm(x(k,:)-C-d(i)*V'); 
    distL1(i)=sum(abs(x(k,:)-C-d(i)*V'));     
end;

 dPQSQ(k) = d(find(distPQSQ == min(distPQSQ)));
 dL2(k) = d(find(distL2 == min(distL2)));
 dL1(k) = d(find(distL1 == min(distL1)));
 
end

plot(U,dL2,'go','LineWidth',3); hold on;
plot(U,dPQSQ,'bo','LineWidth',3); hold on;
plot(U,dL1,'mo','LineWidth',3); hold on;
plot(U,U,'r--','LineWidth',5); hold on;