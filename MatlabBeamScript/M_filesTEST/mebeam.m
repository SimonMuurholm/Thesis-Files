function Me = mebeam(Xe,He)
%.........................................................
% mebeam: 
%   Creates the element mass matrix Me of a beam 
%   element.
%
% Syntax:
%   Me = mebeam(Xe,me)
%
% Input:
%   Xe   : initial coordinates   = [x1  y1
%                                   x2  y2].
%   He  : element properties [E A I (G As) m]
%
% Output:
%   Me   : element mass matrix.
%
% Version 1.0    11.02.13
% Version 2.0    12.03.13
%.........................................................

% Form initial element (column) vector a0.
a0 = (Xe(2,:)-Xe(1,:))';  % element vector
L = sqrt(a0'*a0);         % element length
n = a0/L;                 % unit vector

me = He(size(He,2));

% Element mass matrix.
Me = me*L/420*...
    [ 140    0      0       70    0      0      ;   
        0  156     22*L      0   54    -13*L    ;
        0   22*L    4*L^2    0   13*L   -3*L^2  ;
       70    0      0      140    0      0      ;
        0   54     13*L      0  156    -22*L    ;
        0  -13*L   -3*L^2    0  -22*L    4*L^2 ];
        
Ae = [n(1) n(2)  0    0    0   0  ;
     -n(2) n(1)  0    0    0   0  ;
       0    0    1    0    0   0  ;
       0    0    0   n(1) n(2) 0  ;
       0    0    0  -n(2) n(1) 0  ;
       0    0    0    0    0   1] ;

Me = Ae'*Me*Ae;
