function res = testSetGenerator( n, m, q, p, mu, f, N)
%testSetGenerator generates test set for comparison of different version of
%PCA. Resulting matrix res is returned as output argument and is written to
%disk with file name created by specified ruler
%   Arguments:
%   n is number of observations.
%   m is number of dimensions.
%   q is number of uniformly distributed dimensions (“true” dimension of space).
%   p is number of contaminated attributes.
%   mu is magnitude of outliers.
%   f is fraction of outliers.
%   N number of dataset

    %Generate uniformly distributed values U(-0.5,0.5)
    res = rand(n,m)-0.5;
    
    %Transform first q columns to U(-10,10)
    res(:, 1:q) = res(:, 1:q)*20;

    %Transform all other columns to Laplacian noise L(0,0.1)
    b = 0.1 / sqrt(2);
    res(:,q+1:m) = - b * sign(res(:,q+1:m)).* log(1- 2* abs(res(:,q+1:m)));
    
    %Select an outliers
    if p>0 && f>0
        %Number of outliers
        k = round(n*f);
        %Get number of cases which become outliers
        ind = randsample(n,k);
        %Change the magnitude of outliers to L(mu,0.1).
        res(ind,q+1:q+p) = res(ind,q+1:q+p) + mu;
    end
    
    %Test set is ready. Form file name
    fname = ['n', num2str(n,'%04.0f'), 'm', num2str(m,'%04.0f'), 'q',...
        num2str(q), 'p', num2str(p), 'mu', num2str(mu,'%02.0f'), 'f'...
        '#', num2str(N,'%03.0f')];
    %Save test set into file
    dlmwrite(fname,res,'\t');
end

