function [inside, newZPoly]=Algorithm5(ZPoly, arcs, ZPoint)
% Function to calculate the point in region problem for complex shapes
% delimited by linear, circular or 3rd order Bezier segments.
% Points are represented in complex notation for easier manipulations
%
% The parameters are:
% ZPoly - matrix of columns of 3 points (complex numbers) to designate segments
%   Start point
%   Center point (for arcs) or first control point for Bezier
%   Second control point for Bezier
%     Center point and Second point can be NaN to distinguish types of arcs

% arcs - vector expressing the sense of arcs
%    0 - no arc for the next point (straight line)
%   +1 - circular arc to next point in the positive sense (CCW)
%   -1 - circular arc to next point in the negative sense (CW)

% ZPoint - array of 1 or more testing points in the complex representation

% inside - array with logical values for the insideness of ZPoints
% newZPoly - the actual set of vertices used (where absent values were passed)
%
% The algorithm allows some random generation of parameters if they are
% absent in the data passed in ZPoly:
%    Center point of arcs if absent
%    Control points in Bezier (not really random but some auto-creation)

% This functions uses 3 other functions for better code organization:
% v=NRAintLin(z1,z2, plotg)
%     For linear segments
% v=NRAintArc(z1,z2,z0,sense,plotg)
%     For circular arcs
% v=NRAintBez(z,plotg)
%     For Bezier segments

% (c) V. Santos, vitor@ua.pt, May 2019, July 2024

ZPoly(1,end+1)=ZPoly(1,1);  %close polyline for easier calculation
arcs(end+1)=nan;

% Add random centers for circular arcs where they do not exist
for n=1:numel(arcs)-1
    if( abs(arcs(n))== 1 && isnan(ZPoly(2,n)) )
        z1=ZPoly(1,n);     % start
        z2=ZPoly(1,n+1);   % end   
        zA=(z1+z2)/2;      % mean point
        zB=(z2-z1)/2;      % vector of chord
        zB=zB/norm(zB);    % versor
        zD=zB*exp(1j*pi/2);%direction of line with center
        rr=(2-4*rand); if abs(rr) < 0.1; rr=rr*5; end
        z0=zA+rr*zD;       % place of center with some "radius" 
        ZPoly(2,n)=z0;     % Adjust center in ZPoly
    end
end

% Add "random" control bezier points where they do not exist
for n=1:numel(arcs)-1
    if( arcs(n)== 2 && ( isnan(ZPoly(2,n)) || isnan(ZPoly(3,n)) )  )
        z1=ZPoly(1,n);       % start
        z2=ZPoly(1,n+1);     % end   
        zA=(z1+z2)/4-2j;     % aux calc
        zB=(z1+z2)/2*1.25+2j;% aux calc
        ZPoly(2,n)=zA;       % Adjust CNTRL1 in ZPoly
        ZPoly(3,n)=zB;       % Adjust CNTRL2 in ZPoly
    end
end

%ZPoly

inside=false(size(ZPoint));

for m=1:numel(ZPoint)           % trying all points
    ZPolyTemp=ZPoly-ZPoint(m);  % temp polygon for a given point
    ZPolyTemp(isnan(ZPolyTemp))=NaN+1i*NaN; %to erase some possible imaginary valid in reals NaN :-( Maybe not needed!
    v=0;                        % initial crossing counter
    for n=1:size(ZPolyTemp,2)-1
        z1=ZPolyTemp(1,n);      % start
        z2=ZPolyTemp(1,n+1);    % end
        z0=ZPolyTemp(2,n);      % circular center or 1st bezier control
        z3=ZPolyTemp(3,n);      % 2nd bezier control  
        switch abs(arcs(n))
            case 0              % linear segment
                v=v+NRAintLin(z1,z2,0);
            case 1              % circular segment
                v=v+NRAintArc(z1,z2,z0, arcs(n),0);
            case 2              % bezier segment
                z=[z1 z0 z3 z2];
                v=v+NRAintBez(z,0);
        end
    end
    if v
        inside(m)=1;
    end
end

newZPoly=ZPoly(:,1:end-1);  %return the final points used (suppressing the last -- repeated one)
