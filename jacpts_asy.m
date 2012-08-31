function [x w] = jacpts_asy(n,a,b)

if a == 0 && b == 0, [x w] = legpts_asy(n); return; end  % call legpts_asy

if n <= 20
    [xbdy wbdy] = jacpts_asy2_bdy(n,a,b,ceil(n/2));  
    [xbdy2 wbdy2] = jacpts_asy2_bdy(n,b,a,floor(n/2));  
    x = [-xbdy2(end:-1:1) ; xbdy];
    w = [wbdy2(end:-1:1)  wbdy];
    return
end

[x w] = jacpts_asy1(n,a,b);

% changeover
nbdy = 10;
bdyidx1 = n-(nbdy-1):n;
bdyidx2 = nbdy:-1:1;

% boundary
[xbdy wbdy] = jacpts_asy2_bdy(n,a,b,nbdy);  
x(bdyidx1) = xbdy;  w(bdyidx1) = wbdy;
if a ~= b
    [xbdy wbdy] = jacpts_asy2_bdy(n,b,a,nbdy);  
end
x(bdyidx2) = -xbdy; w(bdyidx2) = wbdy;