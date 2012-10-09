function [x w] = legpts_asy(n)
% Computes the legendre nodes and weights using asymptotic formulae and
% Newton's method. 

% interior
[x w] = legpts_asy1(n);   

% changeover
nbdy = 10;
bdyidx1 = n-(nbdy-1):n;
bdyidx2 = nbdy:-1:1;

% boundary
[xbdy wbdy] = legpts_asy2_bdy(n,nbdy);  
x(bdyidx1) = xbdy;  w(bdyidx1) = wbdy;
x(bdyidx2) = -xbdy; w(bdyidx2) = wbdy;
