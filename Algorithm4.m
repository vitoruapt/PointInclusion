function M=Algorithm4(P,A)
% P - polygon vertices
% A - list of points to test
% M - vector of logicals with true for inclusion
%
% Optimized version where no divisions are performed :-))

% To speed up in case of many points we first remove
% all testing points that exceed limits of vertices

inc=find(~(A(1,:)>max(P(1,:)) | A(1,:)<min(P(1,:)) | A(2,:)>max(P(2,:)) | A(2,:)<min(P(2,:))) );
M=false(1,size(A,2));
P=[P P(:,1)]; %add an extra vertex similar to the first for easier indexing
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
        if Y1==0 || Y2==0 %this rarely occurs (so can skip to speed up?)
            Ip=Ip*0.5;
        end
        IC=IC+Ip;
    end
    
    if IC ~= 0 
        M(n)=true;
    end
end
end