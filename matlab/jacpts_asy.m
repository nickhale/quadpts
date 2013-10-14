function [x, w] = jacpts_asy(n,a,b)

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

if ( a == 0 && b == 0)
    % Call legpts_asy:
    [x, w] = legpts_asy(n);
    return
end  

% TODO: Special cases, a, b = \pm 0.5?

if ( n <= 20 )
    % Use only the boundary formula:
    [xbdy, wbdy] = jacpts_asy2_bdy(n, a, b, ceil(n/2));  
    [xbdy2, wbdy2] = jacpts_asy2_bdy(n, b, a, floor(n/2));  
    x = [-xbdy2(end:-1:1) ; xbdy];
    w = [wbdy2(end:-1:1),  wbdy];
    return
end

% Call the main (interior) routine:
[x, w] = jacpts_asy1(n,a,b);

% Changeover:
nbdy = 10;
bdyidx1 = n-(nbdy-1):n;
bdyidx2 = nbdy:-1:1;

% Call the boundary routine:
[xbdy, wbdy] = jacpts_asy2_bdy(n, a, b, nbdy);  
x(bdyidx1) = xbdy;  
w(bdyidx1) = wbdy;

% Call for the negative points (if a ~= b):
if ( a ~= b )
    [xbdy, wbdy] = jacpts_asy2_bdy(n, b, a, nbdy);  
end
x(bdyidx2) = -xbdy; 
w(bdyidx2) = wbdy;

end