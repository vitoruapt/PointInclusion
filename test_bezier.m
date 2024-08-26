% Used to test bezier and generate some images for paper
% Mar, April 2019, July 2024
% (c) vitor@ua.pt, 2019,2024

close all
clear

syms P1x P2x P3x P4x P1y P2y P3y P4y real 
P=[P1x P2x P3x P4x
   P1y P2y P3y P4y];

B=bezier(P)

Q=[0  2   1   5  5
   0  5  -4   5  2
   ];

Q=[0  2   1   -1  3
   0  3  -4   0  2
   ];
Q=circshift(Q,1);


Q=[0  2   1   3  3
   0  3  -14  3  2
   ];


%a symmetric "S" curve
D=0; %control points deviation from extremes
Offx=0; %Horizontal offset
Offy=0; %Vertical offset
fEnd=1;

%   start   CTRL1     CTRL2       end
Q=[ Offx    D+Offx    Offx+fEnd-D    fEnd+Offx  
    Offy    Offy+fEnd    Offy        fEnd+Offy ];

%Q=circshift(Q,1); % comment for other orientation
%Q(1,:)=2*Q(1,:);

a1=[3/4;1/2];
a2=[5/4;1/4];
a1=[0.9;0.6];
a2=[1.1;0.4];
a3=[1;0.5];
a4=[0.5;0.2];


Q=Q-a3;

% z0=0; z1=2+2j; z2=-1j; z3=2+1j;
% Q=[ 0 2 0 2
%     0 2 -1 1];
% Q=Q-a1;


P1x=Q(1,1);
P2x=Q(1,2);
P3x=Q(1,3);
P4x=Q(1,4);
P1y=Q(2,1);
P2y=Q(2,2);
P3y=Q(2,3);
P4y=Q(2,4);


Bb=collect(subs(B));
DBb=collect(diff(subs(B),'t'));
DDBb=collect(diff(diff(subs(B),'t'),'t')); %second derivative if you need it.
% %-------Intermediate--------------------
% latex(Bb)
% latex(DBb)
% latex(sym(Q)')
% return
% %-----------------------

cc=coeffs(Bb(2));  % coefficients of polynomial y(t)
cc=flip(cc); %coeffs gives them in reverse order :-(
tc=roots(cc); %numeric way to calculate roots in Matlab
tc=double(tc);  %convert from symbolic to double

%idx=find(isreal(tc));  %this may fails for very small imaginary parts...
idx=find(abs(imag(tc))<1e-16);  %...very small imaginary parts are discardable
tcr=tc(idx);  %array with only the real ones
tcr=real(tcr); %eliminate remains of neglegible imaginary parts

%calculate the values of derivatives at the intersections
mydel=1e-5; %a value to operate near the left neighborhood and resolve ambiguity due to null derivative

t=tc-mydel;
valB=double(subs(Bb(1)))';
valDB=double(subs(DBb(2)))';

%plot the polynomial curve
t=linspace(0,1,500);
G=double(subs(B));
plot(G(1,:),G(2,:),'b-');
hold on
grid on
plot(0,0,'og');


%plot the intersections
t=tcr(:)'-mydel;
GI=double(subs(B));
plot(GI(1,:),GI(2,:),'r*');
%plot with squares the valid NRA (negative real axis)
idxP=find(GI(1,:)<0);  %or <=?
plot(GI(1,idxP),GI(2,idxP),'bs'); axis equal

%label all intersection NRA and other
lbl=repmat(['t='],size(tcr,1),1);
txt=[lbl num2str(tcr,4)]; 
text(GI(1,:)+0.02,GI(2,:)-0.05,txt);


%adjust some axis for figure export
%axis([-0.9 1.1 -0.6 0.4]) %a1
%axis([-1.1 0.9 -0.4 0.6]) %a2

%
%valDB=vpa(subs(DBb(2),5));

%t=t(idxP);

%plot( Q(1,:), Q(2,:), 'k*')

% % x0=fEnd+Offx;
% % y0=Offy;
% % r=fEnd;
% % t=linspace(pi/2, pi, 200);
% % Gx=x0+r*cos(t);
% % Gy=y0+r*sin(t);
% % plot(Gx,Gy,'b-'); axis equal
IC=-(GI(1,:)'< 0).*sign(valDB')   %shouldn't this be <= and not only < ???
disp('    t       x(t)   dy/dt    IC');
%vpa([tcr G(1,:)'  valDB'  IC],4)

vpa([tc valB(1,:)'  valDB'  IC],4)


return
%------------

Bb=simplify(B)
DBb=simplify(diff(B,'t'))
