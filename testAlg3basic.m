% Script for simple vairous teste on Algorithm 3 (of supported Paper)
close
clear

% a CCW polygon
Pccw=[
     -2     2     2    -2
     -2    -2     2     2
];

% the same but CW polygon
Pcw=[
     -2    -2     2     2
     -2     2     2    -2
];

% A polygon with a segment over NRA
PoverNRA = [
     -2     2     2    -2
     -2    -2     0     0
];

% A rotated square
Plosang= [
    -2  -1  0  -1
     0  -1  0   1
     ];

% A trap geometry :-)
Ptrap= [
    -2  -1  0  0  -1  -1
     0  -1  0  2   0   1
     ];


P=Pcw;
P=Ptrap;

A=[-2 ; 0];  %on left border
A=[2 ; 0];  %on right border
A=P; %all the vertices
A=[-0.5  -1.5
    1.5   0.5];

M=Algorithm3(P,A);

%Plot results
h=fill(P(1,:),P(2,:),'y');
h.FaceColor=0.5*[1 1 1];
h.FaceAlpha=0.2;
hold on, grid on, axis square, axis equal%, axis(0.5*[-5 5 -5 5])
plot(A(1,M),A(2,M),'g.','MarkerSize',20); 
plot(A(1,~M),A(2,~M),'r.','MarkerSize',20);

h=line( [P(1,6) 0.5],[P(2,6) 2.5]);
h.Color=[0 0 0 ];
h.LineStyle=":";
