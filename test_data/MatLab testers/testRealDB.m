function res = testRealDB( data, fun, numOfInt, opt, delta, maxPC )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    res = zeros(maxPC, 4);
    N = size(data,1)*size(data,2);
    %Loop of #PC
    for pc = 1:maxPC
        res(pc,1) = pc;
        
        %Estimate time for one calculation
        tic;
        [V,U,C] = pcaPQSQ(data, pc, fun, 'optimize', opt, ...
            'numofintervals', numOfInt, 'intshrinkage', delta);
        %Calculate projections
        rest = bsxfun(@plus,U*V',C);
        tim = toc;
        
        if tim<0.0001
            nRep = 10000;
        elseif tim<1
            nRep = round(1/tim);
        else
            nRep = 1;
        end
        
        if nRep>1
            tic;
            for p = 1:nRep
                [V,U,C] = pcaPQSQ(data, pc, fun, 'optimize', opt, ...
                    'numofintervals', numOfInt, 'intshrinkage', delta);
                %Calculate projections
                rest = bsxfun(@plus,U*V',C);
            end
            tim = toc;
        end
        
        res(pc,2) = tim/nRep;
        res(pc,3) = sum(abs(data(:)-rest(:)))/N;
        res(pc,4) = var(data(:)-rest(:));
    end
end

