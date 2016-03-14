function [V,U,C] = pcaL1(x, ncomp, varargin)

verbose=0;
useJavaImplementation = 0;


for i=1:length(varargin)
    if strcmp(varargin{i},'verbose')
        verbose = varargin{i+1};
    end
    if strcmp(varargin{i},'javacode')
        useJavaImplementation = 1;
    end
end

if ~useJavaImplementation
    varargin{length(varargin)+1} = 'potential';
    varargin{length(varargin)+1} = @L1;
    [V,U,C] = pcaPQSQ(x,ncomp,varargin{:});
else
    varargin{length(varargin)+1} = 'potential';
    varargin{length(varargin)+1} = 'L1';
    [V,U,C] = pcaPQSQ(x,ncomp,varargin{:});
end


end

