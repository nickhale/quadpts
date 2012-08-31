function [x w v t] = jacpts_asy2_bdy(n,a,b,npts)

% if npts > ceil((n+1)/2)
%     error('NPTS must be <= N/2');
% end

% Useful constants
rho = n + .5*(a + b + 1); rho2 = n + .5*(a + b - 1);

smallk = min(30,npts);
% Use GLR for finding the first bessel roots
jk = besselroots(a,min(npts,smallk));
% use asy formula for larger ones (See NIST 10.21.19, Olver 74 p247)
if npts > smallk
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
if any(t > pi/2), warning('jacpts_bdy:theta','Theta > pi/2'); end

% [tB1 A2] = asy2_higherterms(a,b);
[tB1 A2 tB2 A3] = asy2_higherterms(a,b);
% [tB1 A2 tB2 A3 tB3 A4] = asy2_higherterms(0,0);

dt = inf; j = 0;
% Newton iteration
while norm(dt,inf) > sqrt(eps)/200
    [vals ders] = feval_asy(n,t,0);    % Evaluate via asymptotic formula
    dt = vals./ders;                   % Newton update
    t = t + dt;                        % Next iterate
    j = j + 1; if j > 10, dt = 0; end  % Bail
end

% ab = a + b; x = cos(t);
% P = .5*(a-b+(ab+2)*x);  Pm1 = 1; 
% Pp = .5*(ab+2);         Ppm1 = 0; 
% for k = 1:n-1
%     A = 2*(k+1)*(k+ab+1)*(2*k+ab);
%     B = (2*k+ab+1)*(a^2-b^2);
%     C = prod(2*k+ab+(0:2)');
%     D = 2*(k+a)*(k+b)*(2*k+ab+2);
%     
%     Pa1 = ( (B+C*x).*P - D*Pm1 ) / A;
%     Ppa1 = ( (B+C*x).*Pp + C*P - D*Ppm1 ) / A;
%     
%     Pm1 = P; P = Pa1;  
%     Ppm1 =  Pp; Pp = Ppa1;
% end
% t = acos(x-P./Pp);

[ignored ders] = feval_asy(n,t,1);     % Evaluate via asymptotic formula

% Constant for weights
if a && b
    if n > 50
        M = min(20,n-1); CW = 1; phi = -a*b/n;
        for m = 1:M
            CW = CW + phi;
            phi = -(m+a)*(m+b)/(m+1)/(n-m)*phi;
            if abs(phi/CW) < eps/100, break, end
        end
    else
        CW = gamma(n+a+1)*gamma(n+b+1)/gamma(n+a+b+1)/factorial(n);
    end
    CW = 2^(a+b+1)*CW;
else
    CW = 2^(a+b+1);
end

% flip
t = t(npts:-1:1); ders = ders(npts:-1:1);% vals = vals(npts:-1:1);
% Revert to x-space
x = cos(t);      w = [CW./ders.^2].';   v = sin(t)./ders;

    function [vals ders] = feval_asy(n,t,flag)
        
        % Useful constants
        A = (.25-a^2);       B = (.25-b^2);
        
        % Evaluate the Bessel functions
        Ja = besselmx(double('J'),a,rho*t,0);
        Jb = besselmx(double('J'),a+1,rho*t,0);
        Jbb = besselmx(double('J'),a+1,rho2*t,0);
        if ~flag
            Jab = besselmx(double('J'),a,rho2*t,0);
        else
            % In the final step, perform accurate evaluation
            Jab = besseltaylor(-t,rho*t);
%             tmpa = Ja;
%             tmpb = Jab;
%             Ja = double(besselj(a,vpa(rho)*t));
%             Jab = double(besselj(a,vpa(rho2)*t));
%             semilogy(cos(t),abs(Ja-tmpa)./sin(t),'b'); hold on
%             semilogy(cos(t),abs(Jab-tmpb)./sin(t),'r'); 
%             Ja = tmpa; Jab = tmpb;
%             title(['a = ' num2str(a) ', b = ', num2str(b)])
%             figure
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
        [nu JK] = meshgrid(-kmax:kmax, z);
        Bjk = besselmx(double('J'),a+nu,JK,0);
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

% j(1) = vpa('2.5574510185965304905903315638721353345562818990684132',64);
% j(2) = vpa('5.6756963202731099020336037475172203722555749045325732',64);
% j(3) = vpa('8.8099925220020443572462957670605946054738411139407751',64);
% j(4) = vpa('11.948059752924419967228033420138143006090687505955351',64);
% jk(1) = 2.5574510185965304;
% jkb(1) = 3.408444291286342e-17;
% jk(2) = 5.6756963202731099;
% jkb(2) = 2.390986187939089e-16;
% jk(3) = 8.8099925220020443;
% jkb(3) = -1.665415490426429e-16;
% jk(4) = 11.9480597529244199672280;
% jkb(4) = -3.841773544987923e-16;


%     function [Ja Jab] = accurateBesselJ(a, tt)
%         % Compute h = rho*t - jk accurately
%         tta = double(single(tt));    ttb = tt - tta;
%         rhoa = double(single(rho));  rhob = rho - rhoa;
%         hi = rhoa*tta;               lo = rhoa*ttb + rhob*tt;
%         h = (hi-jk) + lo;            h2 = h - tt;
%         kmax = min(ceil(abs(log(eps)/log(norm(h2,inf)))),30);
%         H = bsxfun(@power,h,0:kmax).';
%         H2 = bsxfun(@power,h2,0:kmax).';
% 
%         % Compute coeffs in Taylor expansions about jk (See NIST 10.6.7)
%         [nu JK] = meshgrid(-kmax:kmax, jk);
%         Bjk = besselmx(double('J'),a+nu,JK,0);        
% %         Bjk(1:min(7,npts),kmax+2) = double(besselj(a+1,vpa(jk(1:min(7,npts)))));
% %         Bjk(:,kmax) = -Bjk(:,kmax+2);
%         nck = abs(pascal(floor(1.25*kmax),1)); nck(1,:) = []; % nchoosek
%         AA = [Bjk(:,kmax+1) zeros(npts,kmax)];  
%         fact = 1;
%         for k = 1:kmax
%             sgn = -1;
%             for l = 0:k
%                 sgn = -sgn;
%                 AA(:,k+1) = AA(:,k+1) + sgn*nck(k,l+1)*Bjk(:,kmax+2*l-k+1);
%             end
%             fact = k*fact;
%             AA(:,k+1) = AA(:,k+1)/2^k/fact;
%         end
% 
%         % Evaluate Taylor series
%         Ja = zeros(npts,1); Jab = zeros(npts,1);
%         for k = 1:npts
%             Ja(k,1) = AA(k,:)*H(:,k);
%             Jab(k,1) = AA(k,:)*H2(:,k);
%         end
%         
%     end