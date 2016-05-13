function res = testRealDBL2( data, maxPC )
    res = zeros(maxPC, 4);
    N = size(data,1)*size(data,2);
    %Loop of #PC
    for pc = 1:maxPC
        res(pc,1) = pc;
        
        %Estimate time for one calculation
        tic;
        c = sum(data)/size(data,1);
        dat = bsxfun(@minus,data,c);
        %Calculate PCs
        [U, S, V] = svds(dat,pc);
        %Calculate projections
        rest = bsxfun(@plus,U*S*V',c);
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
                c = sum(data)/size(data,1);
                dat = bsxfun(@minus,data,c);
                %Calculate PCs
                [U, S, V] = svds(dat,pc);
                %Calculate projections
                rest = bsxfun(@plus,U*S*V',c);
            end
            tim = toc;
        end
        
        res(pc,2) = tim/nRep;
        res(pc,3) = sum(abs(data(:)-rest(:)))/N;
        res(pc,4) = var(data(:)-rest(:));
    end
end

