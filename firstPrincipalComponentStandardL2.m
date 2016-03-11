function [V,U,C,explainedUariance,MDS] = firstPrincipalComponentStandardL2(x, varargin)

verbose=0;
eps=0.001;


for i=1:length(varargin)
    if strcmp(varargin{i},'verbose')
        verbose = varargin{i+1};
    end
    if strcmp(varargin{i},'eps')
        eps = varargin{i+1};
    end
end

C = mean(x);
C1 = C;
V = zeros(1,size(x,2));
V1 = zeros(1,size(x,2));
U = zeros(size(x,1),1);
U1 = U;
N = size(x,1);
mn = mean(x);
stdev = std(x);

Xc = zeros(size(x,1),size(x,2));
for i=1:size(x,1) 
    Xc(i,:)= x(i,:)-mn; 
end;

XU = zeros(1,size(x,2));
SU = 0;
SU2 = 0;

for k=1:size(x,2)
    V(k) = rand()*stdev(k);
end

count = 1;

V = V/norm(V);

while(count<1000)

U1 = U;
C1 = C;
    
for i=1:size(x,1)
    U(i)=sum(Xc(i,:).*V)/sum(V.*V);
    %st=st+U(i);
    SU2=SU2+U(i)*U(i);
    %Xt=Xt+x(i,:)*U(i);
    XU=XU+Xc(i,:)*U(i);
end

for k=1:size(x,2)
    %det = st2*N-st*st;
    %V1(k) = (Xt(k)*N-mn(k)*N*st)/det;
    %A(k) = (mn(k)*N*st2-st*Xt(k))/det;
    V1(k) = XU(k)/SU2;
end

if(V1(1)<0) V1=-V1;
end
V1 = V1/norm(V1);

delta = (norm(V-V1))/norm(V1);

if(delta<eps)
    V = V1;
    break;
end

V=V1;


if verbose 
    mds = 0;
    deltaU = (norm(U-U1))/norm(U1);
    deltaC = (norm(C-C1))/norm(C1);
    for i=1:size(x,1)
        diff = x(i,:)-U(i)*V;
        mds=mds+sum(diff.*diff);
    end
    mds = mds/N;
    MDS(count) = mds;
    %display(sprintf('%i: Delta=%f, DeltaU=%f, DeltaC=%f, MDS=%f', count, delta, deltaU, deltaC, mds));
    display(sprintf('%i: DeltaV=%f, DeltaU=%f, MDS=%f', count, delta, deltaU, mds));
end

count = count+1;

end

explainedUariance = var(U);

if verbose
    display(sprintf('Fraction of explained variance: %f', explainedUariance/sum(var(x))));
end
