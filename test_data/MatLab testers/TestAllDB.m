%preallocation
maxPC = 10;
meths = 11;
res = zeros(maxPC, 4, meths);
%L2
res(:,:,1) = testRealDBL2( t1, maxPC );
names = cell(1,10);
names{1} = 'L2';
k = 2;
res(:,:,k) = testRealDB( t1, @L1, 5, 0, 1, maxPC );
names{k} = 'L1,5,0,1';
k = k+1;
res(:,:,k) = testRealDB( t1, @L1, 5, 1, 1, maxPC );
names{k} = 'L1,5,1,1';
k = k+1;
res(:,:,k) = testRealDB( t1, @L1, 5, 2, 1, maxPC );
names{k} = 'L1,5,2,1';
k = k+1;
res(:,:,k) = testRealDB( t1, @L1, 10, 0, 1, maxPC );
names{k} = 'L1,10,0,1';
k = k+1;
res(:,:,k) = testRealDB( t1, @L1, 10, 1, 1, maxPC );
names{k} = 'L1,10,1,1';
k = k+1;
res(:,:,k) = testRealDB( t1, @L1, 20, 0, 1, maxPC );
names{k} = 'L1,20,0,1';
k = k+1;
res(:,:,k) = testRealDB( t1, @L1, 20, 1, 1, maxPC );
names{k} = 'L1,20,1,1';
k = k+1;
res(:,:,k) = testRealDB( t1, @L1, 2, 0, 1, maxPC );
names{k} = 'L1,2,0,1';
k = k+1;
res(:,:,k) = testRealDB( t1, @L1, 2, 0, 0.75, maxPC );
names{k} = 'L1,2,0,0.75';
k = k+1;
res(:,:,k) = testRealDB( t1, @L1, 2, 0, 0.5, maxPC );
names{k} = 'L1,2,0,0.5';
