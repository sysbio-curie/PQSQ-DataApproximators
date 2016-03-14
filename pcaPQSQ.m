function [V,U,C] = pcaPQSQ(x, ncomp, varargin)

verbose=0;
intervals = defineIntervals(x,5);
useJavaImplementation = 0;

for i=1:length(varargin)
    if strcmp(varargin{i},'verbose')
        verbose = varargin{i+1};
    end
    if strcmp(varargin{i},'intervals')
        intervals = varargin{i+1};
    end
    if strcmp(varargin{i},'javacode')
        useJavaImplementation = 1;
    end
    if strcmp(varargin{i},'potential')
        potential_function_handle = varargin{i+1};
    end
end

if ~useJavaImplementation
xwork = x;
C = PQSQ_Mean(x,intervals,potential_function_handle);

for i=1:ncomp
    if verbose
        display(sprintf('Component %i',i));
    end
    
    [Vi,Ui] = firstPrincipalComponentPQSQ(xwork, intervals, potential_function_handle, 'mean', C, varargin);
    
    V(i,:) = Vi;
    U(:,i) = Ui;
    
    %%%%
    % You might want to use PQSQ-based projections: this should increase precision?
    %%%%%
    for j=1:size(x,1)
        proj = 0; s1 = 0; s2 = 0;
        for k=1:size(x,2)
           s1=s1+(xwork(j,k)-C(k))*Vi(k);
           s2=s2+Vi(k)*Vi(k);
        end
        proj = s1/s2;
        xwork(j,:) = xwork(j,:) - C - proj*Vi;
    end
    
end

else % java implementation
    javaclasspath({'VDAOEngine.jar'});
    pca = vdaoengine.analysis.PCAMethodPQSQ;
    pca.setData(x);
    if strcmp(potential_function_handle,'L1')
        pca.PQSQpotential = vdaoengine.analysis.PQSQPotential.getTrimmedLinearPQSQPotential(x);
    end
    if strcmp(potential_function_handle,'L2')
        pca.PQSQpotential = vdaoengine.analysis.PQSQPotential.getTrimmedQuadraticPQSQPotential(x);
    end
    if strcmp(potential_function_handle,'LSQRT')
        pca.PQSQpotential = vdaoengine.analysis.PQSQPotential.getTrimmedSqrtPQSQPotential(x);
    end
    
    pca.PQSQpotential.numberOfIntervals = size(intervals,2);
    pca.PQSQpotential.intervals = intervals;
    pca.PQSQpotential.computeABCoefficients();
    
    pca.verboseMode = 0;
    pca.calcBasis(ncomp);
    C = pca.getBasis().a0;
    C = C';
    U = pca.pointProjections;
    U = U';
    U = U(:,1:ncomp);
    V = pca.getBasis().basis;
end

end

