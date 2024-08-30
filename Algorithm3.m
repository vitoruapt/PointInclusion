function M=Algorithm3(P,A)
% P - polygon vertices
% A - list of points to test
% M - vector of logicals with true for inclusion
% To speed up in case of many points, first we remove
% all testing points that exceed limits of vertices:

%inc=1:size(A,2);
inc=find(~(A(1,:)>max(P(1,:)) | A(1,:)<min(P(1,:)) | A(2,:)>max(P(2,:)) | A(2,:)<min(P(2,:))) );
M=false(1,size(A,2));
P=[P P(:,1)]; %add an extra vertex similar to the first for easier indexing

for n=inc  %test all points that are worth (inside the bounding box of the polygon)
    IC=0;
    Q=P-repmat(A(:,n),1,size(P,2)); %translate to origin
    for m=1:size(P,2)-1
        Y1=Q(2,m);
        Y2=Q(2,m+1);
        DY=Y1-Y2;
        if ( sign(Y1)==sign(Y2) )
            continue  %does not intersect NRA
        end
        X1=Q(1,m);
        X2=Q(1,m+1);
        DX=X1-X2;
        t=Y1/DY;  %safe because DY is not zero!
        
        tDX=t*DX;
        
        if(tDX==X1) %over the border (consider inside)
            IC=1; %irrelevant, simply for the counter below
            break %interrupt the for loop and avoid next calculations
        end

        if (t>=0 && t <=1 && tDX>X1 ) %|| ( (t==0 || t==1) && t*DX==X1)  % must be >X1 and not >=X1 to exclude over the line 
            Ip=sign(DY);
            if t==0 || t==1 %this rarely occurs (intersection at the end of the segment)
                Ip=Ip*0.5;
            end
            IC=IC+Ip;
        end
    end
    
    if IC ~= 0   %total intersection count must be non null to be included
        M(n)=true;
    end
end
end