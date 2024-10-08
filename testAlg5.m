% Script to test poly region point inclusion - i.e. Algorithm 5 from supporting paper
% Define regions by vertices and type of arcs that connect them
% to form the poly region.
% Random points are generated, tested and represented.
% You can create as many points as desired and try zooming the
% plot to confirm the inclusion
%
% (c) V. Santos, vitor@ua.pt, May 2019, July 2024


% ---- general layout of image
close
axis equal
clear
%line([-4 4], [0 0]); line([ 0 0], [-4 4]);
axis off
hold on
sFact=1; %just a default

%% ======= FIRST PART FOR ISOLATED TESTs OF NRA inclusion --
% SKIP to next Matlab Section for full region tests

% z1=4-8*rand + 1i*(4-8*rand); 
% z2=4-8*rand + 1i*(4-8*rand);
% z0=NaN; %to calculate automatically with a random generator

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

%  z0=0; z1=2+2j; z2=-1j; z3=2+1j;
%  z=[z0 z1 z2 z3];
% % a= 1+j/2; %3/4+j/2; %a1=[3/4;1/2]; a2=[5/4;1/4]; a1=[0.9;0.6]; a2=[1.1;0.4]; a3=[1;0.5]; a4=[0.5;0.2];
% a=1.12+0.5i
% vb=NRAintBez(z-a,1), return ;       %just testing


%% ================ SECOND PART CONCERNED WITH FULL POLYREGIONS ===

%This example generates absolutely random polyregions
% Results may be difficult to understand visually and some may be non-sense!
% Hence, try first the examples below

NN=7;

% ----- random polyregion (generates only linear and circular.
%  To generate Bezier you should give additional contro points...
ZPoly=nan(3,NN);
ZPoly(1,:)=4-8*rand(1,NN) + 1i*(4-8*rand(1,NN)); 
%arcs=(1-2*(rand(1,NN)>0.5)); %generates only arcs with positive and negative senses.
arcs=(1-2*(rand(1,NN)>0.5)); mask=rand(1,NN)<0.6; arcs(mask)=0; %generates mixed arcs and linear ( 60% linear)

arcs(1)=-1;
%------------------------------------------

%---- random testing points % see below a more general solution
% MM=3000;
% ZPoint=10-20*rand(1,MM) + 1i*(10-20*rand(1,MM));
% ZPoint=ZPoint/6;

%another way using uniform spacing
% KK=75;
% xx=linspace(-1.6,1.6,KK); yy=linspace(-1.5,1.5,KK);
% [XX,YY]=meshgrid(xx,yy);
% ZPoint=XX+j*YY; ZPoint=ZPoint(:);

%-------------------------------



% == Now follow several region to test the algorithm
% Simply (un)comment the lines of ZPoly, arcs and sFact in each case

%% Example 1 -- particular case of poly line 1
% ZPoly = [ -3-2j  0.5-j  1+2j
%            NaN    NaN  -1+0j
%            NaN    NaN   NaN
%         ];
% arcs=[0 1 -1];
% sFact=2.5;  %aux scale factor for range of test points


%% Example 2 --- Nice particular case :-)
% ZPoly=[
%   -0.6461 - 3.8155i  -0.5638 + 2.7258i   0.0787 - 0.4238i  -3.6338 - 2.3731i   2.9074 + 3.6398i   1.5969 - 2.7054i   2.2940 + 2.0809i   3.4734 + 0.0556i   0.3233 - 1.6164i   1.3943 + 0.8956i
%   -2.2726 - 0.5239i  -0.5104 + 1.0963i  -1.3980 - 2.1214i  -0.5838 + 0.8734i   3.3613 + 0.2381i   1.2259 - 0.2074i   2.4903 + 0.8392i   1.7462 - 0.4938i   1.1295 - 0.4758i  -0.2822 - 1.1757i
%       NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i
% ];
% arcs=[-1     1     1    -1    -1     1    -1    -1     1     1 ];
% sFact=2;  %aux scale factor for range of test points


%% Example 3
%  ZPoly=[
%   -2.3882 + 0.4710i  -2.2797 + 1.8528i   0.0966 + 1.6826i   2.0024 + 3.8855i  -3.5619 - 2.9007i  -0.5562 - 1.0885i  -2.0188 - 1.5156i  -0.6938 - 3.5914i   0.3940 - 2.6554i   1.8645 - 2.4714i
%   -3.8175 + 1.2784i      NaN + 0.0000i  -0.4617 + 4.0915i      NaN + 0.0000i  -1.4880 - 2.9417i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i
%       NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i
% ];
% arcs=[ -1     0     1     0    -1     0     0     0     0     0 ];


%% Example 4
% ZPoly = [ 0      2+1j   0.5-1j  -0-1.5j -2    1+1.5j
%           2+2j    NaN   NaN      NaN   NaN   NaN
%           0-1j    NaN   NaN      NaN   NaN   NaN
%          ];
% arcs=[2 -1 1 -1 2 0];
% sFact=2;  %aux scale factor for range of test points

%% Example 5
 % ZPoly=[
 %   0.0000 + 0.0000i   2.0000 + 1.0000i   0.5000 - 1.0000i   0.0000 - 2.0000i  -2.0000 + 0.0000i   1.0000 + 2.0000i   0.0000 + 0.0000i
 %   2.0000 + 2.0000i   0.2071 + 0.7822i      NaN + 0.0000i   0.1494 + 0.1494i  -0.2500 - 1.5000i      NaN + 0.0000i   0.0000 + 0.0000i
 %   0.0000 - 1.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i  -0.6250 + 3.2500i      NaN + 0.0000i   0.0000 + 0.0000i
 %   ];
 % arcs=[2 1 0 -1 2 0];
 % sFact=1.25;  %aux scale factor for range of test points

%% Example 6

%   ZPoly=[
%    0.0000 + 0.0000i   2.0000 + 1.0000i   0.5000 - 1.0000i   0.0000 - 1.5000i  -2.0000 + 0.0000i   1.0000 + 2.0000i   0.0000 + 0.0000i
%    2.0000 + 2.0000i   0.9662 + 0.2129i      NaN + 0.0000i  -0.2647 + 0.2304i  -0.2500 - 1.5000i      NaN + 0.0000i   0.0000 + 0.0000i
%    0.0000 - 1.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i  -0.6250 + 3.2500i      NaN + 0.0000i   0.0000 + 0.0000i
% ];
% arcs=[2 1 0 -1 2 0];
% sFact=1.25;  %aux scale factor for range of test points


%% Example 7

% ZPoly=[
%    0.0000 + 0.0000i   2.0000 + 1.0000i   0.5000 - 1.0000i   0.0000 - 1.5000i  -2.0000 + 0.0000i   1.0000 + 2.0000i   0.0000 + 0.0000i
%    2.0000 + 2.0000i   0.3040 + 0.7095i  -0.3593 - 0.6407i  -0.5303 - 0.1237i  -0.2500 - 1.5000i      NaN + 0.0000i   0.0000 + 0.0000i
%    0.0000 - 1.0000i      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i  -0.6250 + 3.2500i      NaN + 0.0000i   0.0000 + 0.0000i
% ];
% arcs=[2 1 1 -1 2 0];
% sFact=1.25;  %aux scale factor for range of test points

%% Example 8
% The complex case of paper

ZPoly =[
 -2.25-1.7i   -1.500 - 1.7000i    0-1i  2.0000 + 1.0000i   0.5000 - 1.0000i   0.0000 - 1.5000i  -2.0000 + 0.0000i   1.0000 + 1.5000i  
   NaN         2.0000 + 2.0000i   NaN   0.2596 + 0.7428i  -0.8340 - 0.1660i   0.0404 + 0.6372i  -0.2500 - 1.6250i      NaN + 0.0000i   
   NaN         0.0000 - 1.0000i   NaN      NaN + 0.0000i      NaN + 0.0000i      NaN + 0.0000i  -0.6250 + 2.9375i      NaN + 0.0000i   
];
arcs=[0 2 0 -1 1 -1 2 0];
sFact=1.2;  %aux scale factor for range of test points

%% Example 9
%- The 3 arcs in a triangluar arrangement - Example 1 of paper

% ZPoly = [ 0      1+0j  1.25+j*0.5*sqrt(3)/2   0.75+j*1.5*sqrt(3)/2     0.25+j*1.5*sqrt(3)/2     -0.25+j*0.5*sqrt(3)/2
%           NaN    1.5         NaN                 0.5+j*sqrt(3)                NaN                    -0.5
%           NaN    NaN         NaN                   NaN                        NaN                     NaN
%   ];
% ZPoly=ZPoly-(0.5+j*sqrt(3)/2);
% arcs=[0 1 0 1 0 1];
% sFact=2.1;

%% End of examples - now follow the creation of points and actual tests
%  finally genertate the test points 

MM=10000;
sFact=max(sFact,1.2);             %just a scale factor to generate enough testing points
minR=sFact*min(real(ZPoly(1,:)));
maxR=sFact*max(real(ZPoly(1,:)));
minI=sFact*min(imag(ZPoly(1,:)));
maxI=sFact*max(imag(ZPoly(1,:)));

ZPoint=minR+(maxR-minR)*rand(1,MM) + 1i*(minR+(maxR-minR)*rand(1,MM));

% the real magic occurs here :-)
[inside, newZPoly]=Algorithm5(ZPoly, arcs, ZPoint);

% Now plot it
plot( real(ZPoint(inside)) , imag(ZPoint(inside)) , '.g', 'MarkerSize', 4);
plot( real(ZPoint(~inside)), imag(ZPoint(~inside)), 'xr','MarkerSize',4);

plot(ZPoly(1,1:end-1),'.k','MarkerSize',10);
plot(ZPoly(1,1),'.b','MarkerSize',20);
plot(ZPoly(1,2),'ok','MarkerSize',10);

%plot(ZPoly(1,1:end-1),'.k');
vs_draw_polyregion(newZPoly, arcs);
title("First point is blue dot and second point is black cicumference")

