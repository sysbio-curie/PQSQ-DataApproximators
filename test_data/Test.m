function res = Test()
    dirList = dir;
    %Select test directories only
    D = length(dirList);
    ind = zeros(1,D);
    k=0;
    for d=1:D
        if dirList(d).isdir==1  %only directories
            if strcmp(dirList(d).name(1),'n')  %and starts from n
                k = k+1;
                ind(d) = k;
            end
        end
    end
    %k is number of tests
    for d=1:D
        if ind(d)>0
            str = oneDirTest(dirList(d).name);
            if ind(d)==1
                res = str;
            else
                res = [res, str];
            end
        end
    end
end

function str = oneDirTest( name )
    str.name = name;
    %Read list of files
    fList = dir(name);
    %Calculate number of files
    D = length(fList);
    k = 0;
    for d = 1:D
        if fList(d).isdir==0
            k=k+1;
        end
    end
    %Create arrays for results
    str.time = zeros(k,1);
    str.error = zeros(k,1);
    %Start individual tests
    k = 0;
    for d = 1:D
        if fList(d).isdir==0
            %down load test data
            k=k+1;
            t1 = dlmread([name, '\', fList(d).name],'\t');
            display(fList(d).name);
            nTimes = 1;
            %Evaluate time and calculate projections onto the space of the
            %first principal components.
            tic;
            for p = 1:nTimes
%L2                res = simplePCA(t1);
                res = PQSQ_PCA(t1);
            end
            str.time(k) = toc/nTimes;
            %Calculate the error as sum of absolute values of the 6-10
            %coordinates
            str.error(k) = sum(sum(abs(res(:,6:end))))/size(t1,1);
        end
    end
    
end