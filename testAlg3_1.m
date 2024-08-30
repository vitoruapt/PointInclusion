%test alg 3
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

P=Pcw;

A=[-2 ; 0];  %on left border
A=[2 ; 0];  %on right border
A=P; %all the vertices

M=Algorithm3(P,A);

%Plot results
h=fill(P(1,:),P(2,:),'y');
h.FaceColor=[0.9 0.9 0.9];
h.FaceAlpha=0.2;
hold on, grid on, axis square, axis equal, axis([-5 5 -5 5])
plot(A(1,M),A(2,M),'g*'); 
plot(A(1,~M),A(2,~M),'r*');

