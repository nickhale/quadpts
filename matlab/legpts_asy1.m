function [x, w, v, ders, t] = legpts_asy1(n,mint)
  
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

% Approximate roots via asymptotic formula [Tricomi, 1950].
isOdd = mod(n,2);
k = (n-2+isOdd)/2+1:-1:1; theta = pi*(4*k-1)/(4*n+2);
x = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);
t = acos(x);

% Minimum value of theta to be compute accurately by interior formula
if ( nargin == 1 )
    mint = t(end-9);                      % By default, ignore last 10
elseif ( mint == 0 )
    mint = t(end-1);                      % Ignore the last point
end
idx = max(find(t<mint,1)-1, 1);           % Index of minimum theta

dt = inf; j = 0;
while norm(dt,inf) > sqrt(eps)/1000       % Newton iteration   
    [vals, ders] = feval_asy(n,t,mint,0); % Evaluate via asymptotic formulae
    dt = vals./ders;                      % Newton update
    t = t - dt;                           % Next iterate
    j = j + 1;
    dt = dt(1:idx-1);
    if ( j > 10 ), dt = 0; end
end
[vals,ders] = feval_asy(n, t, mint, 1);  % Once more for luck
t = t - vals./ders;                      % Newton update
    
% Revert to x-space and compute weights
x = cos(t);
w = 2./ders.^2;
v = sin(t)./ders;

% Flip using symetry for negative nodes
if ( isOdd )
    x = [-x(end:-1:2) x].';
    w = [w(end:-1:2) w];
    v = [v(end:-1:2) v].';
    ders = [ders(end:-1:2) ders].';
else
    x = [-x(end:-1:1) x].';
    w = [w(end:-1:1) w];
    v = [v(end:-1:1) v].';
    ders = [ders(end:-1:1) ders].';
end

end


function [vals, ders] = feval_asy(n,t,mint,flag)
M = 20;  % Maximum number of terms

% Coefficients in expansion.
c = cumprod((1:2:2*M-1)./(2:2:2*M));
d = cumprod((1:2:2*M-1)./(2*n+3:2:2*(n+M)+1));
c = [1 c.*d];  
% How many actually required?
R = (8/pi)*c./(2*sin(mint)).^(.5:M+1)/10;
R = R(abs(R)>eps); M = length(R); c = c(1:M);

% Constant out the front
ds = -1/8/n; s = ds; j = 1;
while abs(ds/s) > eps/100
    j = j+1;
    ds = -.5*(j-1)/(j+1)/n*ds;
    s = s + ds;
end
p2 = exp(s)*sqrt(4/(n+.5)/pi);
% Stirling's coefficients for Gamma function
g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
     5246819/75246796800 -534703531/902961561600 ...
     -4483131259/86684309913600 432261921612371/514904800886784000];
fn = sum(g.*[1 cumprod(ones(1,9)./n)]);
fn5 = sum(g.*[1 cumprod(ones(1,9)./(n+.5))]);
C = p2*(fn/fn5);

% Some often used vectors/matrices
onesT = ones(1,length(t));
onesM = ones(M,1);
M05 = transpose((0:M-1)+.5);
onesMcotT = onesM*cot(t);
M05onesT = M05*onesT;

% The main terms
alpha = onesM*(n*t) + M05onesT.*(onesM*(t-.5*pi));
cosAlpha = cos(alpha);
sinAlpha = sin(alpha);

% Accurate computation of cos(alpha(1,:)) via Taylor series
if ( flag )
    k = numel(t):-1:1;
    rho = n+.5;
    % Accurate computation of rho*theta-(k-.25)*pi
    ta = double(single(t));    tb = t - ta;
    hi = rho*ta;               lo = rho*tb;
    pia = double(single(pi));
    pib = -8.742278000372485e-08;   % pib = pi - pia;
    dh = (hi-(k-.25)*pia)+lo-(k-.25)*pib; 
    % Initialise
    tmp = 0; sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh;
    for j = 0:20                    % Sum terms in Taylor expansion
        dc = sgn*DH/fact;
        tmp = tmp + dc;
        sgn = -sgn;
        fact = fact*(2*j+3)*(2*j+2);
        DH = DH.*dh2;
        if norm(dc,inf) < eps/1000, break, end
    end
    tmp(2:2:end) = -tmp(2:2:end);  % Get the sign right
    tmp = sign(cosAlpha(1,2)*tmp(2))*tmp;
    cosAlpha(1,:) = tmp;           % Update cosAlpha
end

% Compute the function evaluation
twoSinT = onesM*(2*sin(t));
denom = cumprod(twoSinT)./sqrt(twoSinT);
vals = C*(c*(cosAlpha./denom));    % Sum up all the terms

% Compute the derivatives
numer = M05onesT.*(cosAlpha.*onesMcotT + sinAlpha) + n*sinAlpha;
ders = -C*(c*(numer./denom));      % Sum up all the terms (dP/dt)

end

