function [tB1 A2 tB2 A3 tB3 A4] = asy2_higherterms(a,b,theta)

% The constants a = alpha and b = beta
A = (.25-a^2); B = (.25-b^2); % These are more useful

% For now, just work on half of the domain
c = pi/2;

N = 30;
N1 = N-1;

% 2nd-kind Chebyshev points and barycentric weights
t = .5*c*(sin(pi*(-N1:2:N1)/(2*N1)).'+1);        
v = [.5 ; ones(N1,1)]; v(2:2:end) = -1; v(end) = .5*v(end);

% The g's
g = A*(cot(t/2)-2./t)-B*tan(t/2);
gp = A*(2./t.^2-.5*csc(t/2).^2)-.5*(.25-b^2)*sec(t/2).^2;
gpp = A*(-4./t.^3+.25*sin(t).*csc(t/2).^4)-4*B*sin(t/2).^4.*csc(t).^3;
g(1) = 0; gp(1) = -A/6-.5*B; gpp(1) = 0;

% B0
B0 = .25*g./t;
B0p = .25*(gp./t-g./t.^2);
B0(1) = .25*(-A/6-.5*B);
B0p(1) = 0;

% A1
A10 = a*(A+3*B)/24;
A1 = .125*gp - (1+2*a)/2*B0 - g.^2/32 - A10;
A1p = .125*gpp - (1+2*a)/2*B0p - gp.*g/16;
A1p_t = A1p./t;
A1p_t(1) = -A/720-A^2/576-A*B/96-B^2/64-B/48+a*(A/720+B/48);

% Make f accurately
fcos = B./(2*cos(t/2)).^2;
f = -A*(1/12+t.^2/240+t.^4/6048+t.^6/172800);
idx = t>1e-3;
ti = t(idx);
f(idx) = A.*(1./ti.^2 - 1./(2*sin(ti/2)).^2);
f = f - fcos;

% Integrals for B1
C = cumsummat(N)*(.5*c); 
D = diffmat(N)*(2/c);
I = (C*A1p_t);
J = (C*(f.*A1));

% B1
tB1 = -.5*A1p - (.5+a)*I + .5*J;
tB1(1) = 0;
B1 = tB1./t;
B1(1) = A/720+A^2/576+A*B/96+B^2/64+B/48+a*(A^2/576+B^2/64+A*B/96)-a^2*(A/720+B/48);

% A2
K = C*(f.*tB1);
A2 = .5*(D*tB1) - (.5+a)*B1 - .5*K;
A2 = A2 - A2(1);

%[tB1b A2b] = asy2_highertermsc(.1,-.3);

if nargout < 3
    if nargin == 3
        % Evaluate for output
        tB1 = bary(theta,tB1,t,v);
        A2 = bary(theta,A2,t,v);
    else
        % Make function for output
        tB1 = @(theta) bary(theta,tB1,t,v);
        A2 = @(theta) bary(theta,A2,t,v);
    end
    return
end

% A2p
A2p = D*A2;
A2p = A2p - A2p(1);
A2p_t = A2p./t;
% Extrapolate point at t = 0
w = pi/2-t(2:end);
w(2:2:end) = -w(2:2:end);
A2p_t(1) = sum(w.*A2p_t(2:end))/sum(w);

% B2
tB2 = -.5*A2p - (.5+a)*(C*A2p_t) + .5*C*(f.*A2);
B2 = tB2./t;
% Extrapolate point at t = 0
B2(1) = sum(w.*B2(2:end))/sum(w);

% A3
K = C*(f.*tB2);
A3 = .5*(D*tB2) - (.5+a)*B2 - .5*K;
A3 = A3 - A3(1);

if nargout < 6
    if nargin == 3
        % Evaluate for output
        tB1 = bary(theta,tB1,t,v);
        A2 = bary(theta,A2,t,v);
        tB2 = bary(theta,tB2,t,v);
        A3 = bary(theta,A3,t,v);
    else
        % Make function for output
        tB1 = @(theta) bary(theta,tB1,t,v);
        A2 = @(theta) bary(theta,A2,t,v);
        tB2 = @(theta) bary(theta,tB2,t,v);
        A3 = @(theta) bary(theta,A3,t,v);
    end
    return
end

% A2p
A3p = D*A3;
A3p = A3p - A3p(1);
A3p_t = A3p./t;
% Extrapolate point at t = 0
w = pi/2-t(2:end);
w(2:2:end) = -w(2:2:end);
A3p_t(1) = sum(w.*A3p_t(2:end))/sum(w);

% B2
tB3 = -.5*A3p - (.5+a)*(C*A3p_t) + .5*C*(f.*A3);
B3 = tB3./t;
% Extrapolate point at t = 0
B3(1) = sum(w.*B3(2:end))/sum(w);

% A3
K = C*(f.*tB3);
A4 = .5*(D*tB3) - (.5+a)*B3 - .5*K;
A4 = A4 - A4(1);

if nargin == 3
    % Evaluate for output
    tB1 = bary(theta,tB1,t,v);
    A2 = bary(theta,A2,t,v);
    tB2 = bary(theta,tB2,t,v);
    A3 = bary(theta,A3,t,v);
    tB3 = bary(theta,tB3,t,v);
    A4 = bary(theta,A4,t,v);
else
    % Make function for output
    tB1 = @(theta) bary(theta,tB1,t,v);
    A2 = @(theta) bary(theta,A2,t,v);
    tB2 = @(theta) bary(theta,tB2,t,v);
    A3 = @(theta) bary(theta,A3,t,v);
    tB3 = @(theta) bary(theta,tB3,t,v);
    A4 = @(theta) bary(theta,A4,t,v);    
end

end