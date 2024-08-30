function v=NRAintBez(z,plotg)
% return NRA intersection contribution for a 3rd order bezier segment
% z expected to have 4 points
%   z0 - starting point
%   z1 - first control point
%   z2 - second control point
%   z3 - ending point
% procedure:
%  1 - find the Roots of B3(t)
%  2 - discard imaginary and those outside interval [0,+1]
%  3 - discard those that generate intersection on the right side of origin ( x>0 )
%  4 - obtain sign of gradient of curve on those intersection points.
%  5 - negate that value: that is the intersection
%  6 - if by any chance t=0 or t=1 half the contibution (touching NRA)

% Still to check: cases of tangency (double roots), especially if they coincide with the end of the line!

% (c) vitor@ua.pt, July 2024

if nargin < 2
    plotg=0;
end

if numel(z) < 4
    z(2:4)=0; %not really useful but at least something is there
end

z0=z(1); z1=z(2); z2=z(3); z3=z(4);

P0=z0; P1=z1; P2=z2; P3=z3; %Using P or z to help....

if plotg
    P = [P0 ; P1; P2; P3] ;
    plot([P0; P3],'*k')  %the anchor points
    plot([P1; P2],'*y')  %the control points
    tt=linspace(0,1,100);  %100 points to draw a bezier
    ZZ=(3*P1-P0-3*P2+P3)*tt.^3 + 3*(P0-2*P1+P2)*tt.^2+3*(P1-P0)*tt + P0;
    plot(ZZ,'g');
    nn=1:4; lbl={};
    for n=nn; lbl{n}=sprintf('z_%d ', n-1); end
    text(real(P)+0.15, imag(P)+0.15, lbl);   
end


% define polygonal coeficients of Bézier
c3=3*P1-P0-3*P2+P3;
c2=3*(P0-2*P1+P2);
c1=3*(P1-P0);
c0=P0;
tc=roots(imag([c3 c2 c1 c0])); %numeric way to calculate roots of y in Matlab

v=0;
for n=1:numel(tc)
    t=tc(n);
    if abs(imag(t)) > 1e-10; continue; end        %discard complex roots
    if ( real(t)<0 || real(t)>1), continue; end   %discard out of range t
    
    xB3=real(c3*t^3+c2*t^2+c1*t+c0);
    if xB3 > 0; continue; end %  xx cordinate is positive (does not cross NRA)
    
    %now we are ready to test the gradient for the sign

    dd=imag(3*c3*t*t + 2*c2*t + c1);
    IC=-sign(dd);
    
    if(t==0 || t==1), v=v/2; end
    
    if plotg
        z= c3*t^3+c2*t^2+c1*t+c0;
        plot(real(z),imag(z),'db');
        plot([P0 P2 ; P1 P3],'--','Color',[1 0.5 0],'LineWidth',1);
    end
    v=v+IC;
end

 