function [x, w, v, t] = jacpts_asy2_bdy(n, a, b, npts)

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

if ( npts > ceil((n+1)/2) )
    error('NPTS must be <= N/2');
end

% Useful constants
rho = n + .5*(a + b + 1); 
rho2 = n + .5*(a + b - 1);

smallk = min(30,npts);
% Use GLR for finding the first bessel roots
jk = besselroots(a,min(npts,smallk));
% use asy formula for larger ones (See NIST 10.21.19, Olver 74 p247)
if ( npts > smallk )
    mu = 4*a^2;
    a8 = 8*((length(jk)+1:npts).'+.5*a-.25)*pi;
    jk2 = .125*a8-(mu-1)./a8 - 4*(mu-1)*(7*mu-31)/3./a8.^3 - ...
          32*(mu-1)*(83*mu.^2-983*mu+3779)/15./a8.^5 - ...
          64*(mu-1)*(6949*mu^3-153855*mu^2+1585743*mu-6277237)/105./a8.^7;
    jk = [jk ; jk2];
end
jk = real(jk(1:npts));

% Approximate roots via asymptotic formula (see Olver 1974)
phik = jk/rho;
t = phik + ((a^2-.25)*(1-phik.*cot(phik))./(8*phik) - ...
    .25*(a^2-b^2)*tan(.5*phik))/rho^2;

% Only first half, x > 0
if ( any(t > 1.1*pi/2) ), warning('jacpts_bdy:theta', 'Theta > pi/2'); end

[tB1, A2, tB2, A3] = asy2_higherterms(a, b);

dt = inf; j = 0;
% Newton iteration
while norm(dt,inf) > sqrt(eps)/200
    [vals, ders] = feval_asy(n,t,0);    % Evaluate via asymptotic formula
    dt = vals./ders;                   % Newton update
    t = t + dt;                        % Next iterate
    j = j + 1; if j > 10, dt = 0; end  % Bail
end

[ignored, ders] = feval_asy(n,t,1);     % Evaluate via asymptotic formula

% Constant for weights
if ( a && b )
    if n > 50
        M = min(20,n-1); CW = 1; phi = -a*b/n;
        for m = 1:M
            CW = CW + phi;
            phi = -(m+a)*(m+b)/(m+1)/(n-m)*phi;
            if ( abs(phi/CW) < eps/100), break, end
        end
    else
        CW = gamma(n+a+1)*gamma(n+b+1)/gamma(n+a+b+1)/factorial(n);
    end
    CW = 2^(a+b+1)*CW;
else
    CW = 2^(a+b+1);
end

% flip
t = t(npts:-1:1); 
ders = ders(npts:-1:1);
% vals = vals(npts:-1:1);

% Revert to x-space
x = cos(t);      
w = (CW./ders.^2).';   
v = sin(t)./ders;

    function [vals, ders] = feval_asy(n,t,flag)
        
        % Useful constants
        A = (.25-a^2);       B = (.25-b^2);
        
        % Evaluate the Bessel functions
        Ja = besselj(a,rho*t,0);
        Jb = besselj(a+1,rho*t,0);
        Jbb = besselj(a+1,rho2*t,0);
        if ~flag
            Jab = besselj(a,rho2*t,0);
        else
            % In the final step, perform accurate evaluation
            Jab = besseltaylor(-t, rho*t);
        end

        % Evaluate functions for recurrsive definition of coefficients.
        gt = A*(cot(t/2)-(2./t)) - B*tan(t/2);
        gtdx = A*(2./t.^2-.5*csc(t/2).^2) - .5*B*sec(t/2).^2;
        tB0 = .25*gt;
        A10 = a*(A+3*B)/24;
        A1 = gtdx/8 - (1+2*a)/8*gt./t - gt.^2/32 - A10;
        tB1t = tB1(t); A2t = A2(t); % Higher terms

        % VALS:
        vals = Ja + Jb.*tB0/rho + Ja.*A1/rho^2 + Jb.*tB1t/rho^3 + Ja.*A2t/rho^4;
        % DERS:
        vals2 = Jab + Jbb.*tB0/rho2 + Jab.*A1/rho2^2 + Jbb.*tB1t/rho2^3 + Jab.*A2t/rho2^4;
        
        % Higher terms (not needed for n > 1000).
        tB2t = tB2(t); A3t = A3(t);
        vals = vals + Jb.*tB2t/rho^5 + Ja.*A3t/rho^6;
        vals2 = vals2 + Jbb.*tB2t/rho2^5 + Jab.*A3t/rho2^6;
        
%         % Super higher terms (not needed for n > 150)
%         tB3t = tB3(t); A4t = A4(t);
%         vals = vals + Jb.*tB3t/rho^7 + Ja.*A4t/rho^8;
%         vals2 = vals2 + Jbb.*tB3t/rho2^7 + Jab.*A4t/rho2^8;
        
        % Constant out the front (Computed accurately!)
        ds = .5*(a^2)/n;
        s = ds; jj = 1;
        while abs(ds/s) > eps/10
            jj = jj+1;
            ds = -(jj-1)/(jj+1)/n*(ds*a);
            s = s + ds;
        end
        p2 = exp(s)*sqrt((n+a)/n)*(n/rho)^a;
        g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
             5246819/75246796800 -534703531/902961561600 ...
             -4483131259/86684309913600 432261921612371/514904800886784000];
        f = @(z) sum(g.*[1 cumprod(ones(1,9)./z)]);
        C = p2*(f(n+a)/f(n))/sqrt(2);

        valstmp = C*vals;
        denom = sin(t/2).^(a+.5).*cos(t/2).^(b+.5);
        vals = sqrt(t).*valstmp./denom;

        % Relation for derivative
        C2 = C*n/(n+a)*(rho/rho2)^a;
        ders = (n*(a-b-(2*n+a+b)*cos(t)).*valstmp + 2*(n+a)*(n+b)*C2*vals2)/(2*n+a+b);
        ders = ders.*(sqrt(t)./(denom.*sin(t)));
        
    end

    function Ja = besseltaylor(t,z)
        kmax = min(ceil(abs(log(eps)/log(norm(t,inf)))),30);
        H = bsxfun(@power,t,0:kmax).';
        % Compute coeffs in Taylor expansions about z (See NIST 10.6.7)
        [nu, JK] = meshgrid(-kmax:kmax, z);
        Bjk = besselj(a+nu,JK,0);
        nck = abs(pascal(floor(1.25*kmax),1)); nck(1,:) = []; % nchoosek
        AA = [Bjk(:,kmax+1) zeros(npts,kmax)];
        fact = 1;
        for k = 1:kmax
            sgn = 1;
            for l = 0:k
                AA(:,k+1) = AA(:,k+1) + sgn*nck(k,l+1)*Bjk(:,kmax+2*l-k+1);
                sgn = -sgn;
            end
            fact = k*fact;
            AA(:,k+1) = AA(:,k+1)/2^k/fact;
        end
        % Evaluate Taylor series
        Ja = zeros(npts,1);
        for k = 1:npts
            Ja(k,1) = AA(k,:)*H(:,k);
        end
    end

end

function jk = besselroots(nu,m)
    % Find m roots of besselj(nu,x)
    
    jk = zeros(m,1); foo = 3;
    if ( nu == 0 )
        xs = 2.404825557695773;
    elseif ( nu > 0 )
        % See Hethcote 1970
        xs = nu + 1.8557*nu^(1/3);
    else
        nu1 = nu + 1;
        % See Piessens 1984
        xs = 2*sqrt(nu+1)*(1 + nu1/4 - 7*nu1^2/96  + 49*nu1^3/1152 - 8363*nu1/276480);
        foo = min(max(2*ceil(abs(log10(nu1))),3),m);
    end

    % The first root
    jk(1) = BesselNewton(nu,xs);
    if ( m == 1 ), return, end
    % The second root
    jk(2) = BesselNewton(nu,jk(1)+.9*pi);
    if (m == 2), return, end
    % Some more roots
    for k = 3:foo
        jk(k) = BesselNewton(nu,jk(k-1)+.99*pi);
    end
    % The rest
    for k = foo+1:m
        jk(k) = BesselNewton(nu,jk(k-1)+pi);
    end

end

function jk = BesselNewton(nu,jk)
    % Use Newton iterations to find roots

    dx = inf; j = 0;
    while ( dx > sqrt(eps)/1000 )
        u = besselj(nu,jk,0);
        du = besselj(nu-1,jk,0)-nu/jk*u;
        dx = u./du;
        jk = jk - dx;
        j = j+1;
        if ( j > 20 ), break; end
    end
    
end
