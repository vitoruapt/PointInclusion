function vs_draw_polyregion(ZPoly, arcstype)
% function to draw a closed polyline made of linear, circular and Bezier arcs
% Zpoly - array of vertices
%   z1 z2 z3 ... zM (z1) -> created if not present
%   c1 c2 c3 ... cM      -> control param 1 - the circle centre z0, if NaN, z0 is calculated randomly; first Bezier control point
%   s1 s2 s3 ... sM      -> control param 2 - second Bezier control point
% arcstype - vector of the type of arc segment (coarly, the number of additional params to define the segment )
%   0     - linear
%   +/-1  - circular with CCW/CW senses
%   2     - 3rd order bezier
%   other values -  not yet considered


if nargin < 2
    arcstype = zeros( size(Zpoly,2));
end

if ZPoly(1,end) ~= ZPoly(1,1)
    ZPoly(1,end+1) = ZPoly(1,1);
    ZPoly(2:end,end+1) = NaN;
end

lColor='k';

for k=1:size(ZPoly,2)-1
    z1=ZPoly(1,k);   % start
    z2=ZPoly(1,k+1); % end
    z0=ZPoly(2,k);   % center/ bezier control point 1
    z3=ZPoly(3,k);   % bezier control point 2
    
    switch(abs(arcstype(k)))
        case 0  %linear
            z = [z1 ; z2] ;
            plot(real(z),imag(z), lColor);
           
        case 1 %circular (must obtain sense as well)
            sense=arcstype(k);
            if isnan(z0)  %must calculate z0 since it was not given
                zA=(z1+z2)/2; %mean point
                zB=(z2-z1)/2; % vector of chord
                zB=zB/norm(zB); %versor
                zD=zB*exp(1j*pi/2); %direction of line with center
                rr=(4-8*rand); if abs(rr) < 0.1; rr=rr*5; end
                z0=zA+rr*zD;  %place of center with some "radius" 
            end
            th1=angle(z1-z0);
            th2=angle(z2-z0);

            r=norm(z1-z0);
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
            tt=linspace(0,1,500);  %500 points to draw an arc
            ZZ = z0 + r*exp(1i*(th1+tt*(th2-th1)));
            plot(ZZ,lColor);

        case 2   %bezier of 3rd order
            tt=linspace(0,1,1000);  %1000 points to draw a bezier
            P0=z1; P1=z0; P2=z3; P3=z2;
            ZZ=(3*P1-P0-3*P2+P3)*tt.^3 + 3*(P0-2*P1+P2)*tt.^2+3*(P1-P0)*tt + P0;
            plot(ZZ,lColor);
    end
                      
            
end            

