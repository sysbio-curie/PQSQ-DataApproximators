function res = PQSQ_PCA( data )
%simplePCA Calculate mean point and 5 PCs. Then calculate projection of
%data points onto space of the first 5 PCs

    [V, U, C] = pcaPQSQ_fast(data, 5, 'optimize', 2, 'potential', @L1);
    %Calculate projections
    res = bsxfun(@plus,U*V',C);
end

