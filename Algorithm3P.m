function M=Algorithm3P(P,A)
% P - polygon vertices
% A - list of points to test
% M - vector of logicals with true for inclusion

% Attempt to paralelize directly Algorithm3
% Division by zero generates either Inf or NaN
% But it turned out not as efficient (perhaps because of the Inf and NaN)
% 

inc=find( ~(A(1,:)>max(P(1,:)) | A(1,:)<min(P(1,:)) | A(2,:)>max(P(2,:)) | A(2,:)<min(P(2,:))) );
M=false(1,size(A,2)); %start with all outside

P1=[P P(:,1)];  %add an extra vertex similar to the first for easier indexing

for n=inc  %test all points that are worth   (inside the bounding box of the polygon)
    Q=P1-repmat(A(:,n),1,size(P1,2)); %translate to origin
    DD=Q(:,1:end-1)-Q(:,2:end);
    DY=DD(2,:);
    DX=DD(1,:);
    t=Q(2,1:end-1)./DY; %care the division (may give Inf or NaN)
    IList=(t>=0 & t <=1 & t.*DX < Q(1,1:end-1)); %list of segments that intersect NRA
    VList=zeros(size(IList));
    VList(IList)=sign(DY(IList));  
    VList(t==0|t==1)= VList(t==0 | t==1)*0.5;  %still need this
    if sum(VList) ~= 0 
        M(n)=true;
    end
    if any(isnan(t)) %is NaN
        M(n)=true;
    end
end
end