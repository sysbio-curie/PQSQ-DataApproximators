function [V,U,C] = pcaPQSQR(x, ncomp, potential_function_handle, varargin)
%pcaPQSQR calculates the first ncomp principal components for dataset in
%matrix x.
%   Input attributes
%       x is n-by-m matrix with n observations and m attributes.
%       ncomp is a number of principal components to calculate.
%       potential_function_handle is handler of majorant function with
%           syntax "function y = L1(x)". Function must be defined for
%           matrix input x.
%       Other input features are optional and can be specified by 
%       Name, Value pairs:
%       'intervals', intervals pair serves to specify user defined intervals.
%           intervals is matrix each row of which corresponds to one 
%           dimension. In each row the first element must be zero.
%       'number_of_intervals', number_of_intervals pair specifies the
%           number of intervals. 
%       'intshrinkage', delta serves to specify delta which is coefficient
%           for intervals shrinkage (see argument delta in defineIntervals.m).
%       'projecttype', 'PQSQ'/'L2' pair serves to specify type of
%           projections to subtract from data after a principal component
%           finding.
%       'eps', positive number is tolerance level for component calculation.
%   Output arguments
%       V is m-by-ncomp matrix of ncomp columns each of which is the
%           principal component.
%       U is n by ncomp matrix of ncomp columns each of which contains
%           lenths of projection of data point (row) onto corresponding
%           principal component.
%       C is 1-by-m the vector of PQSQ mean (coordinates of central point).
%
%Examples
%  x = rand(100,5); %Generate 100 randomly uniformly distributed points in 
%                   %5D
%  Calculate two principal components for 5 intervals L1 potential function
%  [V,U,C] = pcaPQSQR(x, 2, @L1);
%
%  Draw data in the space of the two first original coordinates.
%  subplot(1,2,1);
%  plot(x(:,1),x(:,2),'ro');
%  Draw data in the space of the first two principal components.
%  subplot(1,2,2);
%  plot(U(:,1),U(:,2),'ro');
%
%  Calcualte the projections of data points onto principal components
%  x1 = bsxfun(@plus,U*V',C);


    intervals = 0;
    number_of_intervals = 5;
    typeOfProjection = 'PQSQ';
    delta = 1;
    eps = 0.001;

    %Process Name, Value pairs
    for i=1:length(varargin)
        if strcmp(varargin{i},'intervals')
            intervals = varargin{i+1};
        elseif strcmp(varargin{i},'number_of_intervals')
            number_of_intervals = varargin{i+1};
        elseif strcmp(varargin{i},'projecttype')
            typeOfProjection = varargin{i+1};
        elseif strcmp(varargin{i},'intshrinkage')
            delta = varargin{i+1};
        elseif strcmp(varargin{i},'eps')
            eps = varargin{i+1};
        end
    end

    %Form potential function structure
    if isscalar(intervals)
        potentialFunction = PQSQR_definePotentialFunction(x,...
            number_of_intervals, potential_function_handle, delta);
    else
        potentialFunction.intervals = [intervals, Inf(size(x,2),1)];
        potentialFunction.sqin = potentialFunction.intervals.^2;
        [potentialFunction.A,potentialFunction.B] = ...
            PQSQR_computeABcoefficients(intervals, potential_function_handle);
    end

    %calculate central point
    C = PQSQR_Mean(x, potentialFunction, eps);

    % initiate xwork like x-C and preallocate matrices for V and U
    xwork = bsxfun(@minus,x,C);
    V = zeros(size(x,2),ncomp);
    U = zeros(size(x,1),ncomp);
        
    for i=1:ncomp
        %Calculate one component
        [Vi, Ui] = onePC_PQSQR(xwork, potentialFunction, varargin{:});
        V(:,i) = Vi;
        U(:,i) = Ui;

        %Subtract projection
        if strcmp(typeOfProjection,'L2')
            xwork = xwork - (xwork*Vi)*Vi';
        else %PQSQ
            xwork = xwork - Ui*Vi';
        end
    end
end
