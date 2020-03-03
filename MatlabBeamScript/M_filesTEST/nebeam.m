function Ne = nebeam(Xe,s)
%.........................................................
% nebeam: 
%   Calculates element displacements in an  
%   elastic beam element.
%
% Syntax:
%   Ne = nebeam(Xe,s)
%
% Input:
%   Xe   : initial coordinates   = [x1    y1
%                                   x2    y2].
%   s    : element coordinate
%
% Output:
%   Ne   : element interpolation matrix.
%
% Date:
%   Version 1.0    13.02.06
%.........................................................

% Calculate element length
a0 = (Xe(2,:)-Xe(1,:))';
L  = sqrt(a0'*a0);

% Define interpolation matrix
Ne = [1-s             0               0 s           0            0
      0   1-3*s^2+2*s^3 (s-2*s^2+s^3)*L 0 3*s^2-2*s^3 (-s^2+s^3)*L];
