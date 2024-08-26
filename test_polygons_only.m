function test_polygons_only()
% Function to test inclusion in polygons only
% It also compares the proposed approach with Matlab own version
% It is prepared to use several testing points (all included ahead in this code):
%  test_inclusion() - the simplest with no optimizations and applied to all points
%  Optimized versions
%       test_inclusion_opt(): no divisions are performed :-))


% (C) vitor@ua.pt, July 2024


close

%fixed sample polygon
P=[
     0     8     8    -2    -2     8     8     2     2     3     3     7     7     0
     0     0    -1    -1    -2    -2    -6    -6     2     2    -5    -5    -3    -3
];
P = P+[-4;4];

%P=20*rand(2,200)-10;  %random polygon
%P=sort(P,2); 
%P=sort(P,1);

A=20*rand(2,100)-10;
A(:,end+1)=[-4 ; 4]; %add a special point for a special polygon :-)
A(:,end+1)=[4 ; 4]; %add a special point for a special polygon :-)
A(:,end+1)=[4 ; 3.5]; %add a special point for a special polygon :-)
A(:,end+1)=[-4 ; 3.5]; %add a special point for a special polygon :-)
A(:,end+1)=[-2 ; 3.5]; %add a special point for a special polygon :-)
A(:,end+1)=[0 ; 2]; %add a special point for a special polygon :-)

% % a special case of many vertices (UK map)
load UK.mat;
P=UK; NP=5000; %axis([-5 2 49 59])
%A=[ (max(P(1,:))-min(P(1,:)))*rand(1,NP)-max(P(1,:)); (max(P(2,:))-min(P(2,:)))*rand(1,NP)-max(P(2,:))];
A=[ -8*rand(1,NP)+3; 49+10*rand(1,NP)];


N=10; %number of times to repeat each test

%  The default -----------

subplot(1,3,1)
fill(P(1,:), P(2,:), 'w');
hold on
%plot(A(1,:),A(2,:),'m.'); 

tic
for n=1:N
    M1=test_inclusion(P,A);
end
t1=toc/N;  %default mean time
n1=nnz(M1);

% The optimised ---------
subplot(1,3,2)
fill(P(1,:), P(2,:), 'w');
hold on
tic
for n=1:N
    Mopt=test_inclusion_opt(P,A);
    %Mopt=test_inclusion_opt_parallel(P,A); %The optimized parallel (works but not more efficient)
    %Mopt=test_inclusion_opt_parallel_v2(P,A); %The optimized parallel v2 (works but not more efficient)
end
top=toc/N;   %optimizes mean time
nop=nnz(Mopt);

% The Matlab implementation -------
subplot(1,3,3)

x = A(1,:)'; y = A(2,:)';
xv= P(1,:)'; yv= P(2,:)';
xv=[xv; xv(1)]; yv=[yv; yv(1)];  %close polygon
tic
for n=1:N
    in = inpolygon(x,y,xv,yv);  %the original matlab  
    %in = inpolygon2(x,y,xv,yv);  %use a modified version form matlabs original
end
t3=toc/N;   %matlab mean time
n3=nnz(in);


%------------ Plot the results -----

% the default with no optimization
subplot(1,3,1)
plot(A(1,M1),A(2,M1),'g.'); 
plot(A(1,~M1),A(2,~M1),'r.'); 
str=sprintf('vs1, t=%f,\nn=%d eff.=%2.1f',t1,n1,t3/t1); title(str)

%The result with optimized version
subplot(1,3,2)
plot(A(1,Mopt),A(2,Mopt),'g.'); 
plot(A(1,~Mopt),A(2,~Mopt),'r.'); 
str=sprintf('vsop, t=%f,\nn=%d eff.=%2.1f',top,nop,t3/top); title(str)

%The matlab version
subplot(1,3,3)
plot(xv,yv,x(in),y(in),'.g',x(~in),y(~in),'.r')
str=sprintf('h&h, t=%f,\nn=%d eff.=%2.1f',t3,n3,t1/t3); title(str)

end

%% 
function M=test_inclusion(P,A)
% P - polygon vertices
% A - list of points to test
% M - vector of logicals with true for inclusion
% To speed up in case of many points we could perhaps remove
% all testing points that exceed limits of vertices
% 
inc=1:size(A,2);
%inc=find(~(A(1,:)>max(P(1,:)) | A(1,:)<min(P(1,:)) | A(2,:)>max(P(2,:)) | A(2,:)<min(P(2,:))) );
M=false(1,size(A,2));
P=[P P(:,1)]; %add a last point to ensure it is closed
for n=inc  %test all points that are worth
    IC=0;
    Q=P-repmat(A(:,n),1,size(P,2)); %translate top origin
    for m=1:size(P,2)-1
        Y1=Q(2,m);
        Y2=Q(2,m+1);
        DY=Y1-Y2;
        if ( sign(Y1)==sign(Y2) )
            continue
        end
        X1=Q(1,m);
        X2=Q(1,m+1);
        DX=X1-X2;
        t=Y1/DY;
        
        tDX=t*DX;
        
        if(tDX==X1) %over the border. 
            IC=1; %irrelevant, simply for the counter below
            break %interrupt the for loop and avoid next calculations
        end

        if (t>=0 && t <=1 && tDX>X1 ) %|| ( (t==0 || t==1) && t*DX==X1)  % must be >X1 and not >=X1 to exclude over the line 
            Ip=sign(DY);
            if t==0 || t==1 %this rarely occurs (so can skip to speed up?)
                Ip=Ip*0.5;
            end
            IC=IC+Ip;
        end
    end
    
    if IC ~= 0 
        M(n)=true;
    end
end
end


%%
function M=test_inclusion2(P,A)
% P - polygon vertices
% A - list of points to test
% M - vector of logicals with true for inclusion
% 
% To speed up in case of many points, first we remove
% all testing points that exceed limits of vertices:
inc=find( ~(A(1,:)>max(P(1,:)) | A(1,:)<min(P(1,:)) | A(2,:)>max(P(2,:)) | A(2,:)<min(P(2,:))) );
M=false(1,size(A,2));

P1=[P P(:,1)];  %add a last point (should check if needed)

for n=inc  %test all points that are worth   
    Q=P1-repmat(A(:,n),1,size(P1,2)); %translate top origin
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

%% 
function M=test_inclusion_opt(P,A)
% P - polygon vertices
% A - list of points to test
% M - vector of logicals with true for inclusion
%
% Optimized version where no divisions are performed :-))

% To speed up in case of many points we first remove
% all testing points that exceed limits of vertices
%inc=1:size(A,2);
inc=find(~(A(1,:)>max(P(1,:)) | A(1,:)<min(P(1,:)) | A(2,:)>max(P(2,:)) | A(2,:)<min(P(2,:))) );
M=false(1,size(A,2));
P=[P P(:,1)]; %add a last point to ensure it is closed and consider all segments
for n=inc  %test all points that are worth
    IC=0;
    Q=P-repmat(A(:,n),1,size(P,2)); %translate top origin
    for m=1:size(P,2)-1
        Y1=Q(2,m);
        Y2=Q(2,m+1);
        DY=Y1-Y2;
        if ( sign(Y1)==sign(Y2) )
            continue
        end
        
        if(DY>0) %( Y1 > Y2)
            if(Y1<0) %( Y1 <0 || Y2 >0)
                continue
            end
        else
            if(Y2<0)% (Y2<0 || Y1 > 0)
                continue
            end
        end
        X1=Q(1,m);
        X2=Q(1,m+1);
        DX=X1-X2;
        
        if sign(DY) > 0
             if Y1*DX < X1*DY
                 continue
             end
        else
             if Y1*DX > X1*DY
                 continue
             end
        end 
        Ip=sign(DY);
        %if Y1==0 || Y2==0 %this rarely occurs (so can skip to speed up?)
        %    Ip=Ip*0.5;
        %end
        IC=IC+Ip;
    end
    
    if IC ~= 0 
        M(n)=true;
    end
end
end


%%
function M=test_inclusion_opt_parallel(P,A)
% P - polygon vertices
% A - list of points to test
% M - vector of logicals with true for inclusion
%
% Optimized version where no divisions are performed :-))
% Variant to try to paralellise all segments at once

% To speed up in case of many points we first remove
% all testing points that exceed limits of vertices
%inc=1:size(A,2);
inc=find(~(A(1,:)>max(P(1,:)) | A(1,:)<min(P(1,:)) | A(2,:)>max(P(2,:)) | A(2,:)<min(P(2,:))) );
M=false(1,size(A,2));
P=[P P(:,1)]; %add a last point to ensure it is closed and consider all segments
for n=inc  %test all points that are worth
    Q=P-repmat(A(:,n),1,size(P,2)); %translate top origin    
    Y1=Q(2,1:end-1);
    Y2=Q(2,2:end);
    DY=Y1-Y2;
    X1=Q(1,1:end-1);
    X2=Q(1,2:end);
    DX=X1-X2;
    noCrossBySgn=(sign(Y1)==sign(Y2));
    noCrossByYA=(Y1>Y2  & (Y1<0 | Y2>0));
    noCrossByYB=(Y1<=Y2 & (Y2<0 | Y1>0));
    noCrossByXA=(sign(DY) > 0 & Y1.*DX < X1.*DY);
    noCrossByXB=(sign(DY) < 0 & Y1.*DX > X1.*DY);   
    Ip= ~(noCrossBySgn | noCrossByYA | noCrossByYB | noCrossByXA | noCrossByXB);
    Ip=Ip.*sign(DY).*(1-0.5*(Y1==0 | Y2==0));
    IC=sum(Ip);      
    if IC ~= 0 
        M(n)=true;
    end
end
end

%%
function M=test_inclusion_opt_parallel_v2(P,A)
% P - polygon vertices
% A - list of points to test
% M - vector of logicals with true for inclusion

% Optimized version where no divisions are performed :-))
% Variant to try to paralellise all segments at once
% Second parallel variant with more optimization

% To speed up in case of many points, we first remove
% all testing points that exceed limits of vertices
%inc=1:size(A,2);
inc=find(~(A(1,:)>max(P(1,:)) | A(1,:)<min(P(1,:)) | A(2,:)>max(P(2,:)) | A(2,:)<min(P(2,:))) );
M=false(1,size(A,2));
P=[P P(:,1)]; %add a last point to ensure it is closed and consider all segments
for n=inc  %test all points that are worth
    Q=P-repmat(A(:,n),1,size(P,2)); %translate top origin    
    Y1=Q(2,1:end-1);
    Y2=Q(2,2:end);
    DY=Y1-Y2;
    X1=Q(1,1:end-1);
    X2=Q(1,2:end);
    DX=X1-X2;
    Y1DX=Y1.*DX;
    X1DY=X1.*DY;
    noCrossBySgn=(sign(Y1)==sign(Y2));
    noCrossByYA=(DY>0  & Y1<0 );
    noCrossByYB=(DY<=0 & Y2<0 );
    noCrossByXA=(sign(DY) > 0 & Y1DX < X1DY);
    noCrossByXB=(sign(DY) < 0 & Y1DX > X1DY);   
    Ip= ~(noCrossBySgn | noCrossByYA | noCrossByYB | noCrossByXA | noCrossByXB);
    Ip=Ip.*sign(DY).*(1-0.5*(Y1==0 | Y2==0));
    IC=sum(Ip);      
    if IC ~= 0 
        M(n)=true;
    end
end
end

