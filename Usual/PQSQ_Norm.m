function norm = PQSQ_Norm(x, intervals, potential_function_handle)
%Calculate the PQSQ Error for one or more points vector.
%  x is n-by-m matrix with rows which correspond to data points.
%  intervals is m-by-K matrix with K values of thresholds which corresponds
%       to each coordinate. intervals(:,K) MUST be equal to Inf.
%  potFunc is potential function handle.
% Return matrix PQSQ potentials for each coordinate of each data point

    [n,m] = size(x);
    norm = zeros(n,m);
    for k=1:m
        norm(:,k) = PQSQ(x(:,k),intervals(k,:),potential_function_handle);
    end
end