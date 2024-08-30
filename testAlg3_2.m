%test alg 3
close
clear

% % a special case of many vertices (UK map)
load UK.mat;
P=UK; NP=5000;
A=[ -8*rand(1,NP)+3; 49+10*rand(1,NP)];
algs={'Alg3','Alg4','Alg4P','Matlab'};
costs=zeros(1,4);
%apply algorithm NN times and calculate mean execution time
NN=100;
for n=1:4
switch n
    case 1
        tic
        for k=1:NN
        M=Algorithm3(P,A);
        end
        costs(n)=toc/NN;
    case 2
        tic
        for k=1:NN
        M=Algorithm4(P,A);
        end
        costs(n)=toc/NN;
    case 3
        tic
        for k=1:NN
        M=Algorithm4P(P,A);
        end
        costs(n)=toc/NN;
    case 4  %matlab's 
        x = A(1,:)'; y = A(2,:)';
        xv= P(1,:)'; yv= P(2,:)';
        tic
        for k=1:NN
            M = inpolygon(x,y,xv,yv);          
        end
        costs(n)=toc/NN;
end

%Plot results
subplot(1,4,n)
h=fill(P(1,:),P(2,:),'w');
h.FaceColor=[0.9 0.9 0.9];
h.FaceAlpha=0.2;
hold on%, grid on, axis equal %, axis([-5 5 -5 5])
plot(A(1,M),A(2,M),'g*'); 
plot(A(1,~M),A(2,~M),'r*');
end


for n=1:4
    subplot(1,4,n)
    str={algs{n}, sprintf('%0.3f',costs(n)), sprintf('%0.1f X', costs(end)/costs(n))};
    title(str)
end
