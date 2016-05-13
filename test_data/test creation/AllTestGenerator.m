nFiles = 100;
n=1000;
m=10;
q=5;
%p=0; %0-3;
mus = [1, 5, 10, 25];
f=0.1;

%clear data without outliers
%create folder name
dirName = ['n', num2str(n,'%04.0f'), 'm', num2str(m,'%04.0f'), 'q',...
        num2str(q), 'p0mu0f', num2str(f*10,'%02.0f')];
%Create folder for files
mkdir(dirName);
dirName=[dirName, '/'];
for nn = 1:nFiles
    testSetGenerator( n, m, q, 0, 0, f, nn, dirName);
end

%Create contaminated sets
for p=1:3
    for mmu=1:length(mus)
        %create folder name
        dirName = ['n',num2str(n,'%04.0f'), 'm', num2str(m,'%04.0f'), 'q',...
        num2str(q), 'p', num2str(p), 'mu', num2str(mus(mmu),'%02.0f'),...
        'f', num2str(f*10,'%02.0f')];
        %Create folder for files
        mkdir(dirName);
        dirName=[dirName, '/'];
        for nn = 1:nFiles
            testSetGenerator( n, m, q, p, mus(mmu), f, nn, dirName);
        end
    end
end
