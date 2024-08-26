
% Script to test poly region point inclusion
% define regions by vertices and type of arcs that connect them
% V. Santos, vitor@ua.pt, May 2019, July 2024



% ---- general lay out of image
close
axis equal
line([-4 4], [0 0]);
line([ 0 0], [-4 4]);
hold on

%% ======= FISRT PART FOR ISOLATED TEST OF NRA inclusion -- SKIP to next Matlab Section for full region tests
z1=4-8*rand + 1i*(4-8*rand); 
z2=4-8*rand + 1i*(4-8*rand);
z0=NaN; %to calculate automatically with a random generator

% ---- Case of tangency
% z1=-4-j; 
% z2=0-j;
% z0=-2-1i*2.5;

% ---- Other testing examples
% z1=2.6-1.19j;
% z2=-1.85-1.18j;
% z0=0.374-1.611j;

% z1=-1.63+0.46j;
% z2=3.84+1.35j;
% z0=1+1.5j;

% -----------------------------

%vc=NRAintArc(z1,z2,z0, 1,1), return ; %just testing
%vl=NRAintLin(z1,z2,1), return ;       %just testing


%% ================ SECOND PART FOR FULL POLYREGIONS


NN=10;

% ----- random polyregion
ZPoly=nan(3,NN);
ZPoly(1,:)=4-8*rand(1,NN) + 1i*(4-8*rand(1,NN)); 
arcs=(1-2*(rand(1,NN)>0.5)); %generates only arcs with positive and negative senses.
arcs=(1-2*(rand(1,NN)>0.5)); mask=rand(1,NN)<0.6; arcs(mask)=0; %generates mixed arcs and linear ( 60% linear)

%arcs(1)=-1;
%------------------------------------------

%---- random testing points
MM=10000;
ZPoint=10-20*rand(1,MM) + 1i*(10-20*rand(1,MM));
ZPoint=ZPoint/1.5;
%-------------------------------



% -- particular case of poly line 1
% ZPoly = [ -3-2j  0.5-j  1+2j
%            NaN    NaN  -1+0j
%            NaN    NaN   NaN
%         ];
% arcs=[0 1 -1];
% --------------------


% --- Nice particular case :-)
% ZPoly=[
%   -0.6461 - 3.8155i  -0.5638 + 2.7258i   0.0787 - 0.4238i  -3.6338 - 2.3731i   2.9074 + 3.6398i   1.5969 - 2.7054i   2.2940 + 2.0809i   3.4734 + 0.0556i   0.3233 - 1.6164i   1.3943 + 0.8956i
%   -2.2726 - 0.5239i  -0.5104 + 1.0963i  -1.3980 - 2.1214i  -0.5838 + 0.8734i   3.3613 + 0.2381i   1.2259 - 0.2074i   2.4903 + 0.8392i   1.7462 - 0.4938i   1.1295 - 0.4758i  -0.2822 - 1.1757i
%       NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i
% ];
% arcs=[-1     1     1    -1    -1     1    -1    -1     1     1 ];
% ---------------------


% ----- POSSIBLY DO NOT NEED TO CHANGE MUCH BELOW

ZPoly(1,end+1)=ZPoly(1,1);  %close polyline
arcs(end+1)=nan;

% Add random centers for cicular arcs where they do not exist
for n=1:numel(arcs)-1
    if( abs(arcs(n))== 1 && isnan(ZPoly(2,n)) )
        z1=ZPoly(1,n);   % start
        z2=ZPoly(1,n+1); % end   
        zA=(z1+z2)/2; %mean point
        zB=(z2-z1)/2; % vector of chord
        zB=zB/norm(zB); %versor
        zD=zB*exp(1j*pi/2); %direction of line with center
        rr=(2-4*rand); if abs(rr) < 0.1; rr=rr*5; end
        z0=zA+rr*zD;  %place of center with some "radius" 
        ZPoly(2,n)=z0 %Adjust center in ZPoly
    end
end

inside=false(size(ZPoint));

for m=1:MM  %trying all points
    ZPolyTemp=ZPoly-ZPoint(m);  %temp polygon for a given point
    v=0;  %initial crossing counter
    for n=1:size(ZPolyTemp,2)-1
        z1=ZPolyTemp(1,n);   % start
        z2=ZPolyTemp(1,n+1); % end
        z0=ZPolyTemp(2,n);   % cicular center or 1st bezier control
        z3=ZPolyTemp(3,n);   % 2nd bezier control  
        switch abs(arcs(n))
            case 0  %linear segment
                vt=NRAintLin(z1,z2,0);
            case 1   %circular segment
                vt=NRAintArc(z1,z2,z0, arcs(n),0);
            case 2  %bezier segment
        end
        v=v+vt;
    end
    if v
        inside(m)=1;
    end
end

plot( real(ZPoint(inside)) , imag(ZPoint(inside)) , '.r');
plot( real(ZPoint(~inside)), imag(ZPoint(~inside)), '.b');

plot(ZPoly(1,1:end-1),'k*');
vs_draw_polyregion(ZPoly, arcs);
