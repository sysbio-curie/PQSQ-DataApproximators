function [V,U,C] = pcaL1(x, ncomp, varargin)

verbose=0;

for i=1:length(varargin)
    if strcmp(varargin{i},'verbose')
        verbose = varargin{i+1};
    end
end


[V,U,C] = pcaPQSQ(x,ncomp,@L1,varargin);


end

