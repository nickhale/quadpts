function [x, w, v, t, ders] = legpts_asy2_bdy(n, npts)

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
rho = n + .5; rho2 = n - .5;

% Roots pf the Bessel function J_0 (Precomputed in Mathematica)
jk = [2.404825557695773     5.520078110286310    8.653727912911012 ...
     11.791534439014281    14.930917708487785   18.071063967910922 ...
     21.211636629879258    24.352471530749302   27.493479132040254 ... 
     30.634606468431975    33.775820213573568].';
if ( npts > 11 )
    % Esimate the larger Bessel roots (See Branders et al., JCP 1981).
    p = ((length(jk)+1:npts).'-.25)*pi;
    aa = [0.000813005721543268 0  0.0245988241803681 0 ...
          0.131420807470708 0 0.0682894897349453];
    bb = [0.00650404577261471 0 0.200991122197811 0 1.16837242570470 0 1 0];
    jk = [jk ; p + polyval(aa,p) ./ polyval(bb,p)]; 
end
jk = jk(1:npts);

% Approximation for Legendre roots (See Olver 1974)
phik = jk/rho;
t = phik + (phik.*cot(phik)-1)./(8*phik*rho^2);

% [tB1 A2] = asy2_higherterms(0,0);
[tB1, A2, tB2, A3] = asy2_higherterms(0,0);
% [tB1 A2 tB2 A3 tB3 A4] = asy2_higherterms(0,0);
dt = inf; j = 0;
% Newton iteration
while ( norm(dt,inf) > sqrt(eps)/200 )
    [vals, ders] = feval_asy(n,t,0);    % Evaluate via asymptotic formula
    dt = vals./ders;                   % Newton update
    t = t + dt;                        % Next iterate
    j = j + 1; 
    if ( j > 10 ), dt = 0; end         % Bail
end
[ignored, ders] = feval_asy(n,t,1);     % Evaluate via asymptotic formula

% flip
t = t(npts:-1:1); ders = ders(npts:-1:1);% vals = vals(npts:-1:1);
% Revert to x-space
x = cos(t);% - sin(t).*vals./ders; 
w = (2./ders.^2).';    v = sin(t)./ders;

    function [vals, ders] = feval_asy(n,t,flag)
        
        % Evaluate the Bessel functions
        Ja = besselj(0, rho*t, 0);
        Jb = besselj(1, rho*t, 0);
        Jbb = besselj(1, rho2*t, 0);
        if ( ~flag )
            Jab = besselj(0, rho2*t, 0);
        else
            % In the final step, perform accurate evaluation
            Jab = besseltaylor(-t, rho*t);
        end
        
        % Evaluate functions for recurrsive definition of coefficients.
        gt = .5*(cot(t) - 1./t);
        gtdt = .5*(-csc(t).^2 + 1./t.^2);
        tB0 = .25*gt;
        A1 = gtdt/8 - 1/8*gt./t - gt.^2/32;
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
        
        % Relation for derivative
        ders = n*(-cos(t).*vals + vals2)./sin(t);
        
        % Common factors
        denom = sqrt(t./sin(t));
        ders = ders.*denom;
        vals = vals.*denom;
       
    end

    function Ja = besseltaylor(t,z)
        kmax = min(ceil(abs(log(eps)/log(norm(t,inf)))), 30);
        H = bsxfun(@power,t,0:kmax).';
        % Compute coeffs in Taylor expansions about z (See NIST 10.6.7)
        [nu, JK] = meshgrid(-kmax:kmax, z);
        Bjk = besselj(nu, JK, 0);
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