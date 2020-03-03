function Se = sbeam(T,X,H,U,ns)
%.........................................................
% sbeam: 
%   Calculates section forces along the beam elements.
%
% Syntax:
%   Sen = sbeam(T,X,u,n)
%
% Input:
%   T    :  topology matrix for bar elements.
%   X    :  node coordinate matrix. 
%   H    :  beam element property matrix. 
%   U    :  element array of nodal displacements
%   ns   :  number of data points along the element
%
% Output:
%   Se   :  stresses in elements
%
%   Version 1.0    13.02.06
%.........................................................

% initialize element section force array
Se = zeros(size(T,1),ns,3);

% Global stiffness matrix by loop over elements
for j = 1:size(T,1)
    
    % Define element arrays
    Xe = X(T(j,1:2),:);       % element node coordinates
    He = H(T(j,3),:);         % element properties
    
    % element displacement vector
    Ue = U(T(j,1:2),:);
    ue = reshape(Ue',6,1);
    
    % element stiffness
    [Ke,Ae] = kebeam(Xe,He);
    
    % local element section forces
    fe = (Ae*Ke*Ae')*(Ae*ue);
    fe([1 3 5]) = -fe([1 3 5]);     % change of sign
    
    % Linear interpolation of section forces [N Q M]
    for i = 1:ns
        s = (i-1)/(ns-1);
        Se(j,i,1:3) = fe(1:3)*(1-s) + fe(4:6)*s;
    end
    
end

