function [x, w, v, t] = jacpts_asy1(n,a,b,mint)

%**************************************************************************
%   This file is part of QUADPTS.
%   Copyright (c) 2012 Nick Hale and Alex Townsend
%   
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%**************************************************************************

% Approximate roots via asymptotic formula.
K = (2*(n:-1:1)+a-.5)*pi/(2*n+a+b+1);
tt = K + 1/(2*n+a+b+1)^2*((.25-a^2)*cot(.5*K)-(.25-b^2)*tan(.5*K));

%% First half, x > 0
t = tt(tt <= pi/2);
if ( nargin == 3 )
    mint = t(end-9);
elseif ( mint == 0 )
    mint = t(end-1); 
end
idx = 1:max(find(t<mint,1)-1,1);

dt = inf; j = 0;
% Newton iteration
while ( norm(dt,inf) > sqrt(eps)/100 )
    [vals, ders] = feval_asy(n, a, b, t, idx, 0); % Evaluate via asy formulae
    dt = vals./ders;                              % Newton update
    t = t + dt;                                   % Next iterate
    j = j + 1;
    dt = dt(idx);
    if ( j > 10 ), dt = 0; end
end
[vals, ders] = feval_asy(n, a, b, t, idx, 1);     % Once more for luck
t = t + vals./ders;

%%
% Constant
if ( a && b )
    if ( n > 50 )
        M = min(20,n-1); C = 1; phi = -a*b/n;
        for m = 1:M
            C = C + phi;
            phi = -(m+a)*(m+b)/(m+1)/(n-m)*phi;
            if ( abs(phi/C) < eps/100 ), break, end
        end
    else
        C = gamma(n+a+1)*gamma(n+b+1)/gamma(n+a+b+1)/factorial(n);
    end
    C = 2^(a+b+1)*C;
else
    C = 2^(a+b+1);
end

%%
% Store
x = cos(t);
w = C./ders.^2;
v = (sin(t)./ders);
t1 = t;

%% Second half, x < 0
tmp = a; a = b; b = tmp;
t = pi - tt(1:(n-length(x)));
if ( nargin == 3) 
    mint = t(10);
elseif ( mint == tt(2) )
    mint = t(2); 
end
idx = max(find(t>mint,1),1):numel(t);

dt = inf; j = 0;
% Newton iteration
while ( norm(dt,inf) > sqrt(eps)/100 )
    [vals, ders] = feval_asy(n, a, b, t, idx, 0); % Evaluate via asy formulae
    dt = vals./ders;                         % Newton update
    t = t + dt;                              % Next iterate
    j = j + 1;
    dt = dt(idx);
    if j > 10, dt = 0; end
end
[vals, ders] = feval_asy(n, a, b, t, idx, 1); % Once more for luck
t = t + vals./ders;                  % Newton update

x = [-cos(t), x].';
w = [C./ders.^2 w];
v = [sin(t)./ders v].';
t = [t, t1].';

function [vals, ders] = feval_asy(n, a, b, t, idx, flag)
M = 20;

% Some often used vectors/matrices
onesT = ones(1,length(t));
onesM = ones(M,1);
MM = transpose(0:M-1);

% Constant out the front
dsa = .5*(a^2)/n; dsb = .5*(b^2)/n; dsab = .25*(a+b)^2/n;
ds = dsa + dsb - dsab; s = ds; j = 1; dsold = ds;
while (abs(ds/s) + dsold) > eps/10 
    dsold = abs(ds/s);
    j = j+1;
    tmp = -(j-1)/(j+1)/n;
    dsa = tmp*dsa*a;
    dsb = tmp*dsb*b;
    dsab = .5*tmp*dsab*(a+b);
    ds = dsa + dsb - dsab;
    s = s + ds;
end
p2 = exp(s)*sqrt(2*pi)*sqrt((n+a)*(n+b)/(2*n+a+b))/(2*n+a+b+1);
% Stirling's coefficients for Gamma function
g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
     5246819/75246796800 -534703531/902961561600 ...
     -4483131259/86684309913600 432261921612371/514904800886784000];
fa = sum(g.*[1 cumprod(ones(1,9)./(n+a))]);
fb = sum(g.*[1 cumprod(ones(1,9)./(n+b))]);
fab = sum(g.*[1 cumprod(ones(1,9)./(2*n+a+b))]);
C = p2*(fa*fb/fab)*2/pi;
C2 = C*(a+b+2*n).*(a+b+1+2*n)./(4*(a+n).*(b+n));

% The sine and cosine terms
alpha = (.5*(2*n+a+b+1+MM))*onesT .* (onesM*t) - .5*(a+.5)*pi;
cosAlpha = cos(alpha);
sinAlpha = sin(alpha);

% if flag
    if idx(1) == 1
        k = numel(t):-1:1;
    else
        k = 1:numel(t);
    end
    ta = double(single(t));   tb = t - ta;
    hi = n*ta;                lo = n*tb + (a+b+1)*.5*t;
    pia = double(single(pi));
    pib = -8.742278000372485e-08; %pib = pi - pia;
    dh = (hi-(k-.25)*pia)+lo-.5*a*pia-(k-.25+.5*a)*pib;
    tmp = 0;
    sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh;
    for j = 0:20
        dc = sgn*DH/fact;
        tmp = tmp + dc;
        sgn = -sgn;
        fact = fact*(2*j+3)*(2*j+2);
        DH = DH.*dh2;
        if norm(dc,inf) < eps/2000, break, end
    end
    tmp(2:2:end) = -tmp(2:2:end);
    tmp = sign(cosAlpha(1,2)*tmp(2))*tmp;
    cosAlpha(1,:) = tmp;
% end

% sin and cos for the derivatives
sint = sin(t);      cost = cos(t);
sinT = onesM*sint;  cosT = onesM*cost;
sinAlpha2 = sinAlpha.*cosT - cosAlpha.*sinT;
cosAlpha2 = cosAlpha.*cosT + sinAlpha.*sinT;

% Compute the PHI function for the values
j = 0:M-2;
vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*n+a+b+j+2);
P1 = [1  cumprod(vec)];
P1(3:4:end) = -P1(3:4:end);
P1(4:4:end) = -P1(4:4:end);
P2 = eye(M);
for l = 0:M-1
    j = 0:(M-l-2);
    vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*n+a+b+j+l+2);
    P2(l+1+(1:length(j)),l+1) = cumprod(vec);
end
PHI = repmat(P1,M,1).*P2;

% Compute the PHI function for the derivatives
j = 0:M-2;
vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*(n-1)+a+b+j+2);
P1 = [1  cumprod(vec)];
P1(3:4:end) = -P1(3:4:end);
P1(4:4:end) = -P1(4:4:end);
P2 = eye(M);
for l = 0:M-1
    j = 0:(M-l-2);
    vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*(n-1)+a+b+j+l+2);
    P2(l+1+(1:length(j)),l+1) = cumprod(vec);
end
PHI2 = repmat(P1,M,1).*P2;

% Compute each term in a FOR loop
S = 0; S2 = 0; 
SC = [ones(1,length(t)); cumprod(onesM(2:end)*(.5*csc(.5*t)))];
halfsec2t = .5*sec(.5*t);
for m = 0:M-1
    l = 0:2:m;
    phi = PHI(m+1,l+1);
    dS1 = phi*SC(l+1,:).*cosAlpha(m+1,:);
    
    phi2 = PHI2(m+1,l+1);
    dS12 = phi2*SC(l+1,:).*cosAlpha2(m+1,:);

    l = 1:2:m;
    phi = PHI(m+1,l+1);
    dS2 = phi*SC(l+1,:).*sinAlpha(m+1,:);
    
    phi2 = PHI2(m+1,l+1);
    dS22 = phi2*SC(l+1,:).*sinAlpha2(m+1,:);

    if m > 10 && norm(dS1(idx) + dS2(idx),inf) < eps/100, break, end
    
    S = S + dS1 + dS2;
    S2 = S2 + dS12 + dS22;
    
    SC(1:m+1,:) = bsxfun(@times,SC(1:m+1,:),halfsec2t);
end

% Multiply by the constant out the front
vals = C*S; S2 = C2*S2;
% Use relation for derivative
ders = (n*(a-b-(2*n+a+b)*cost).*vals + 2*(n+a)*(n+b)*S2)/(2*n+a+b)./sint;
% Compute the denominator
denom = 1./real(sin(t/2).^(a+.5).*cos(t/2).^(b+.5));
% Compute the values and the derivatives
vals = vals.*denom;
ders = ders.*denom;

end

end