%testing algorithm 3,4 of the supporting paper (and compare with Matlab
%default)
close
clear

% a special case of many vertices (UK map)
load UK.mat;
P=UK; NP=5000;

%random points for the UK map
A=[ -8*rand(1,NP)+3; 49+10*rand(1,NP)];

%Only points very close to the border of the UK. Uncomment to test
%A=P; noise=randn(size(A))/1000; A=A+noise;

%---------------------------------------
algs={'Alg3','Alg4','Alg [8]'};
costs=zeros(1,numel(algs));
countin=zeros(1,numel(algs));
countout=zeros(1,numel(algs));

%apply algorithm NN times and calculate mean execution time
NN=10;  %10 for faster results but used 100 for more sustained benchmark
for n=1:numel(algs)
switch n
    case 1
        tic
        for k=1:NN
        M=Algorithm3(P,A);
        end
        costs(n)=toc/NN;
        countin(n)=nnz(M);
        countout(n)=nnz(~M);
    case 2
        tic
        for k=1:NN
        M=Algorithm4(P,A);
        end
        costs(n)=toc/NN;
        countin(n)=nnz(M);
        countout(n)=nnz(~M);
    case 3  %matlab's 
        x = A(1,:)'; y = A(2,:)';
        xv= P(1,:)'; yv= P(2,:)';
        tic
        for k=1:NN
            M = inpolygon(x,y,xv,yv);          
        end
        costs(n)=toc/NN;
        countin(n)=nnz(M);
        countout(n)=nnz(~M);
end

%Plot results
subplot(1,numel(algs),n)
h=fill(P(1,:),P(2,:),'w');
h.FaceColor=[0.9 0.9 0.9];
h.FaceAlpha=0.2;
hold on, axis tight,axis equal%, grid on, axis equal %, axis([-5 5 -5 5])
plot(A(1,M),A(2,M),'g*','MarkerSize',3); 
plot(A(1,~M),A(2,~M),'r*','MarkerSize',3);
end


for n=1:numel(algs)
    subplot(1,numel(algs),n)
    str={algs{n}, sprintf('%0.3f s',costs(n)), sprintf('%0.1f X', costs(end)/costs(n))};
    title(str)
    axis off
    xlabel(sprintf('%d/%d points in/out',countin(n),countout(n) ))
end
