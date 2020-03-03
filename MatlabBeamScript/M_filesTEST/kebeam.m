function [Ke,Ae] = kebeam(Xe,He)
%.........................................................
% kebeam: 
%   Creates the element stiffness matrix of
%   a beam element.
%
% Syntax:
%   Ke = kebeam(Xe,He)
%   [Ke,Ae] = kebeam(Xe,He)
%
% Input:
%   Xe  : nodal coordinates
%   He  : element properties [E A I (G As) m]
%
% Output:
%   Ke   : element stiffness matrix.
%
% Version 1.0    13.02.06
%.........................................................

% Unit directional vector
a0 = (Xe(2,:)-Xe(1,:))';
L = sqrt(a0'*a0);
n = a0/L;

% Element properties
E = He(1);
A = He(2);
I = He(3);

% Shear flexibility parameter
% if size(He,2)>4 && He(4)>0
%     P = 12*E*I/(He(4)*He(5)*L^2);
% else
%     P = 0;
% end

P = 0;

% Local element stiffness matrix
Ke = [E*A/L                 0                 0 -E*A/L                 0                 0
          0  12*E*I/L^3/(1+P)   6*E*I/L^2/(1+P)      0 -12*E*I/L^3/(1+P)   6*E*I/L^2/(1+P)
          0   6*E*I/L^2/(1+P) (4+P)*E*I/L/(1+P)      0  -6*E*I/L^2/(1+P) (2-P)*E*I/L/(1+P)
     -E*A/L                 0                 0  E*A/L                 0                 0
          0 -12*E*I/L^3/(1+P)  -6*E*I/L^2/(1+P)      0  12*E*I/L^3/(1+P)  -6*E*I/L^2/(1+P)
          0   6*E*I/L^2/(1+P) (2-P)*E*I/L/(1+P)      0  -6*E*I/L^2/(1+P) (4+P)*E*I/L/(1+P) ];
        
% Transformation matrix
Ae = [n(1) n(2)  0    0    0   0
     -n(2) n(1)  0    0    0   0
       0    0    1    0    0   0
       0    0    0   n(1) n(2) 0
       0    0    0  -n(2) n(1) 0
       0    0    0    0    0   1];

% Global element stiffness matrix
Ke = Ae'*Ke*Ae;