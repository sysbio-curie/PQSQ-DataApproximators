function [ output_args ] = showComponentsInPCA(x, Vs, C, alpha)
% This function shows a set of vectors, projected into the PCA space of X

[pc,proj] = pca(x);
mn = mean(x);
Cm = C - mn;

plot(proj(:,1),proj(:,2),'ko'); hold on;

plot([1.2*min(proj(:,1)) 1.2*max(proj(:,1))],[0 0],'k-'); hold on;
plot([0 0],[1.2*min(proj(:,2)) 1.2*max(proj(:,2))],'k-'); hold on;


C1 = sum(pc(:,1).*Cm');
C2 = sum(pc(:,2).*Cm');
plot(C1,C2,'ro','MarkerSize',20); hold on;

for i=1:size(Vs,1)
    p1 = project(C+alpha*Vs(i,:),mn,pc(:,1),pc(:,2));
    p2 = project(C-alpha*Vs(i,:),mn,pc(:,1),pc(:,2));
    plot([p1(1) p2(1)],[p1(2) p2(2)],'r-');
end

axis equal;
end

function y=project(x,mn,pc1,pc2)
    xm = x-mn;
    y(1) = sum(xm.*pc1');
    y(2) = sum(xm.*pc2');
end

