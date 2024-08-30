function v=NRAintLin(z1,z2, plotg)
% return NRA intersection contribution for a linear segment
% z1 - starting point (in complex format)
% z2 - ending point (in complex format)

% This version is not the optimized one for linear segments but the
% difference in performance for few points is not relevant.

% (c) vitor@ua.pt, July 2024

if nargin < 3
    plotg=0;
end

if nargin < 2
    z2=z1;
end

v=0;


if plotg
    P = [z1 ; z2] ;
    plot(real(P),imag(P),'m');
    plot(P,'*k')
    nn=1:2; lbl={};
    for n=nn; lbl{n}=sprintf('z_%d ', n); end
    text(real(P)+0.15, imag(P)+0.15, lbl);   
end

Y1=imag(z1);
Y2=imag(z2);
DY=Y1-Y2;
if ( sign(Y1)==sign(Y2) ), return, end
X1=real(z1);
X2=real(z2);
DX=X1-X2;
t=Y1/DY;
if t>=0 && t <=1 && t*DX>X1  % must be <X1 and not <=X1 to include over the line 
    v=sign(DY);
    if t==0 || t==1 %segment extreme
       v=v*0.5;
    end
    
    if plotg
        z=z1+t*(z2-z1);
        plot(real(z),imag(z),'db');
    end
end

 