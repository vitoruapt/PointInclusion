function B=bezier(P)
%function to calculate a bezier curve based on its control points
%Returns a symbolic value for easier usage
% (c) vitor@ua.pt, July 2024

n=size(P,2)-1;
B=0;
syms t
for i=0:n
    c=factorial(n)/factorial(n-i)/factorial(i);
    B=B+c*(1-t)^(n-i)*t^i*P(:,i+1);
end
    