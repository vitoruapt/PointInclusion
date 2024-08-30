function M=Algorithm4P(P,A)
% P - polygon vertices
% A - list of points to test
% M - vector of logicals with true for inclusion
%
% Optimized version where no divisions are performed :-))
% Variant to try to paralellise all segments at once

inc=find(~(A(1,:)>max(P(1,:)) | A(1,:)<min(P(1,:)) | A(2,:)>max(P(2,:)) | A(2,:)<min(P(2,:))) );
M=false(1,size(A,2));
P=[P P(:,1)]; 
for n=inc
    Q=P-repmat(A(:,n),1,size(P,2)); %translate to origin    
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