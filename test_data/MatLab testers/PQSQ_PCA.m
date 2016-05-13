function res = PQSQ_PCA( data )
%simplePCA Calculate mean point and 5 PCs. Then calculate projection of
%data points onto space of the first 5 PCs

    [V, U, C] = pcaPQSQ(data, 5, @L1, 'optimize', 2, 'potential');
    %Calculate projections
    res = bsxfun(@plus,U*V',C);
end

