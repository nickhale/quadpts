function [x, w] = legpts_asy(n)

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

% Computes the lLegendre nodes and weights using asymptotic formulae and
% Newton's method. 

% Changeover point:
nbdy = 10;
bdyidx1 = n-(nbdy-1):n;
bdyidx2 = nbdy:-1:1;

if ( n <= 2*nbdy )
    % Use only the boundary formula:
    [xbdy, wbdy] = legpts_asy2_bdy(n, ceil(n/2));  
    [xbdy2, wbdy2] = legpts_asy2_bdy(n, floor(n/2));  
    x = [-xbdy2(end:-1:1) ; xbdy];
    w = [wbdy2(end:-1:1),  wbdy];
    return
end

% Interior algorithm:
[x, w] = legpts_asy1(n);  

% Boundary algorithm:
[xbdy, wbdy] = legpts_asy2_bdy(n,nbdy);  

% Output the result:
x(bdyidx1) = xbdy;  
w(bdyidx1) = wbdy;
x(bdyidx2) = -xbdy; 
w(bdyidx2) = wbdy;

end
