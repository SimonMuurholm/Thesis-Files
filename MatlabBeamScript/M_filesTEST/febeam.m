function fe = febeam(p,Xe)
%.........................................................
% febeam: 
%   Creates equivalent element load vector for
%   distributed element load.
%
% Syntax:
%   fe = febeam(p,Xe)
%
% Input:
%   p   :  element load array
%   Xe  :  element nodal initial coordinates
%
% Output:
%   fe  :  equivalent element load vector
%
%.........................................................

% Unit directional vector
a0 = (Xe(2,:)-Xe(1,:))';
L = sqrt(a0'*a0);
n = a0/L;

p = [n(1) n(2)  
    -n(2) n(1)]*p';

% Local equivalent element load vector
fe = [p(1)*L/2  p(2)*L/2  p(2)*L^2/12  p(1)*L/2  p(2)*L/2  -p(2)*L^2/12]';

% Transformation matrix
Ae = [n(1) n(2)  0    0    0   0;
     -n(2) n(1)  0    0    0   0;
       0    0    1    0    0   0;
       0    0    0   n(1) n(2) 0;
       0    0    0  -n(2) n(1) 0;
       0    0    0    0    0   1];

% Global equivalent element load vector
fe = Ae'*fe;