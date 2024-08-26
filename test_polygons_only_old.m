% Early version...
function test_polygons_only_old()
clear
close

%fixed sample polygon
P=[
     0     8     8    -2    -2     8     8     2     2     3     3     7     7     0
     0     0    -1    -1    -2    -2    -6    -6     2     2    -5    -5    -3    -3
];
P = P+[-4;4];

P=20*rand(2,50)-10;  %random polygon
%P=sort(P,2); 
%P=sort(P,1);

subplot(1,2,1)
fill(P(1,:), P(2,:), 'w');

%----------------------------------
hold on
A=20*rand(2,10000)-10;
A(:,end+1)=[-4 ; 4]; %add a special point for a special polygon :-)

tic
M=test_inclusion(P,A);
t1=toc;
n1=nnz(M);


% tic
% M=test_inclusion2(P,A);
% t1b=toc;
% n1b=nnz(M);
% str=sprintf('vs2, t=%f, n=%d',t1b,n1b); title(str)



% PV=fliplr(P);
% tic
% M=test_inclusion(PV,A);
% t2=toc
% subplot(1,3,2)
% plot(A(1,M),A(2,M),'g.'); 
% plot(A(1,~M),A(2,~M),'r.'); 

subplot(1,2,2)

x = A(1,:)'; y = A(2,:)';
xv= P(1,:)'; yv= P(2,:)';
xv=[xv; xv(1)]; yv=[yv; yv(1)];  %close polygon
tic
in = inpolygon(x,y,xv,yv);
t3=toc;
n3=nnz(in);


subplot(1,2,1)
plot(A(1,M),A(2,M),'g.'); 
plot(A(1,~M),A(2,~M),'r.'); 
str=sprintf('vs1, t=%f, n=%d eff.=%2.1f',t1,n1,t3/t1); title(str)

subplot(1,2,2)
plot(xv,yv,x(in),y(in),'.g',x(~in),y(~in),'.r')
str=sprintf('h&h, t=%f, n=%d eff.=%2.1f',t3,n3,t1/t3); title(str)

end

%-----------------------------------------------------------
function M=test_inclusion(P,A)
% P - polygon vertices
% A - list of points to test
% M - vector of logicals with true for inclusion
% To speed up in case of many points we could perhaps remove
% all testing points that exceed limits of vertices
% 
inc=find(~(A(1,:)>max(P(1,:)) | A(1,:)<min(P(1,:)) | A(2,:)>max(P(2,:)) | A(2,:)<min(P(2,:))) );
M=false(1,size(P,2));
P=[P P(:,1)]; %add a last point to ensure it is closed
for n=inc  %test all points that are worth
    IC=0;
    Q=P-repmat(A(:,n),1,size(P,2)); %translate top origin
    for m=1:size(P,2)-1
        Y1=Q(2,m);
        %Y2=Q(2,mod(m,size(Q,2))+1);
        Y2=Q(2,m+1);
        DY=Y1-Y2;
        if ( sign(Y1)==sign(Y2) )
            continue
        end
        X1=Q(1,m);
        X2=Q(1,mod(m,size(Q,2))+1);
        DX=X1-X2;
        t=Y1/DY;
        if t>=0 && t <=1 && t*DX<X1  % must be <X1 and not <=X1 to include over the line 
            Ip=sign(DY);
             %if t==0 || t==1 %this rarely occurs (so can skip to speed up?)
             %    Ip=Ip*0.5;
             %end
            IC=IC+Ip;
        end
    end
    
    if IC ~= 0 
        M(n)=true;
    end
end
end


% =========================================
function M=test_inclusion2(P,A)
% P - polygon vertices
% A - list of points to test
% M - vector of logicals with true for inclusion
% 
% To speed up in case of many points we could perhaps remove
% all testing points that exceed limits of vertices

inc=find( ~(A(1,:)>max(P(1,:)) | A(1,:)<min(P(1,:)) | A(2,:)>max(P(2,:)) | A(2,:)<min(P(2,:))) );
M=false(1,size(P,2));
P1=[P P(:,1)];  %add a last point (should check if needed)

for n=inc  %test all points that are worth   
    Q=P1-repmat(A(:,n),1,size(P1,2)); %translate top origin
    %DD=diff(Q,1,2); %does not waht you want, use next
    DD=Q(:,1:end-1)-Q(:,2:end);
    DY=DD(2,:);
    DX=DD(1,:);
    t=Q(2,1:end-1)./DY;   %care the Inf (?)
    IList=(t>=0 & t <=1 & t.*DX< Q(1,1:end-1)); %list of segments that intersect NRA
    VList=zeros(size(IList));
    VList(IList)=sign(DY(IList));  
    VList(t==0|t==1)= VList(t==0 | t==1)*0.5;  %still need this
    if sum(VList) ~= 0 
        M(n)=true;
    end
end
end