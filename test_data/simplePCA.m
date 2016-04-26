function res = simplePCA( data )
%simplePCA Calculate mean point and 5 PCs. Then calculate projection of
%data points onto space of the first 5 PCs
    %Calculate mean and subtract it
    c = sum(data)/size(data,1);
    dat = bsxfun(@minus,data,c);
    %Calculate PCs
    [U, S, V] = svds(dat,5);
    %Calculate projections
    res = bsxfun(@plus,U*S*V',c);
end

