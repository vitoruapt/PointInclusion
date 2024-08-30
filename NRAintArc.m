function v=NRAintArc(z1,z2,z0,sense,plotg)
% return NRA intersection contribution for a circular arc
% z1 - starting point
% z2 - ending point
% z0 - center point (can be absent or NaN for a random generation - Useful only when
%                    drawing locally otherwise may give a result that does not match
%                    drawing somewhere else since this random is not the same as there!!!!)
% sense - +1 CCW (positive),  -1 CW (negative)
% plotg - if 1, plot graphs
% function can also plot a few things for debug. Comment on real usage

% (c) vitor@ua.pt, July 2024

if nargin < 5
    plotg = 0;
end

if nargin < 4
    sense = 1;
end

%Generate a random center if non existing
if nargin < 3 || isnan(real(z0))
    zA=(z1+z2)/2;    % mean point
    zB=(z2-z1)/2;    % vector of chord
    zB=zB/norm(zB);  % versor
    zD=zB*exp(1j*pi/2); %direction of line with center
    rr=(4-8*rand); if abs(rr) < 0.1; rr=rr*5; end
    z0=zA+rr*zD;     % place of center with some "radius" 
end

th1=angle(z1-z0);
th2=angle(z2-z0);

%thetasinicial=180*[th1 th2]/pi;

r=norm(z1-z0);
y0=imag(z0);
x0=real(z0);

v=0;

if th2 <= th1 && sense == 1
    if th2 < 0
        th2=th2+2*pi;
    else
        th1=th1-2*pi;
    end
else
    if th1 <=th2 && sense == -1
        if th1 < 0 
            th1=th1+2*pi;
        else
            th2=th2-2*pi;
        end
    end
end

Dth=th1-th2;


if plotg
    %----- plot points and draw the arc: the old way
    %q = linspace(th1,th2,100); ZZ = z0 + r*exp(1i*q);

    tt=linspace(0,1,100);
    P=[z0 z1 z2];
    hold on
    plot(P,'*')

    ZZ = z0 + r*exp(1i*(th1+tt*(th2-th1)));
    plot(ZZ);
    nn=1:3;
    lbl={};
    for n=nn; lbl{n}=sprintf('z_%d ', n-1); end
    text(real(P)+0.15, imag(P)+0.15, lbl);

    % z = [z1 ; z0; z2;] ; plot(real(z),imag(z));  %plot lines to points
    % z3=z0 + r*exp(1i*(th1+0.5*(th2-th1))); plot(real(z3),imag(z3),'*k'); % plot central point
    %-----------------------

end

%this could be earlier to accelerate procedure
if abs(y0) >= r; return; end  %preserve the = otherwise bad assessment!

aa=asin(-y0/r);

% Solution of: x = sin(q)
% q = asin(x) +/- 2 k pi , k=0,1,2
% q = (pi - asin(x)) +/- 2 k pi = -asin(x) +/- (2k+1)pi , k=0,1,2
% 5 potential solutionsin the interval +/-2pi
% check only those that result in 0 <= tc <= 1

tc1=(th1-aa)/Dth;
tc2=tc1-2*pi/Dth;
tc3=tc1+2*pi/Dth;
tc4=(th1+aa)/Dth-pi/Dth;
tc5=(th1+aa)/Dth+pi/Dth;

allt=round([tc1 tc2 tc3 tc4 tc5],12);  %12 decimal point precision
validtmask=(allt>=0 & allt<=1);

vt=allt(validtmask);

if plotg
    zz=z0 + r*exp(1i*(th1+vt*(th2-th1)));
    plot(real(zz),imag(zz),'sm');
end

for t=allt
    aT=th1-t*Dth;
    IC=0;   
    if t >=0 && t <= 1
        xx=x0+r*cos(aT);
        if xx < 0        % left side of x axis
            IC=-sign(-Dth*cos(aT));
        end
        
%         if aT> -pi/2 && aT<pi/2
%          IC=sign(Dth);
%         else
%             if aT >=pi/2 && aT<=3*pi/2
%                 IC=-sign(Dth);
%             end
%         end
    end
    if t==0 || t==1
        IC=0.5*IC;
    end
    v=v+IC;
end

%ths=[thetasinicial ; 180*[th1 th2]/pi]

