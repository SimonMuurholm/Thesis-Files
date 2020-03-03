function M = mbeam(M,T,X,H,dof)
%.........................................................
% kbeam: 
%   Creates and assembles geometric stiffness matrix
%   of elastic beam elements.
%
% Syntax:
%   M = kgbeam(M,T,X,G)
%  
% Input:
%   M   :  existing global mass matrix.
%   T   :  topology matrix for beam elements.
%   X   :  nodal coordinate matrix.
%   H   :  element property matrix H = [E A I (G As) m];
%   dof :  degrees of freedom per node.
%
% Output:
%   M   :  new global geometric stiffness matrix.
%
% Version 1.0    11.02.13
%.........................................................

% Global stiffness matrix by loop over elements
for j = 1:size(T,1)
    
    % Define element arrays
    Xe = X(T(j,1:2),:);       % element node coordinates
    He = H(T(j,3),:);         % element properties
    
    % element mass matrix
    Me = mebeam(Xe,He);
    
    % element stiffness into global format
    M = assem(M,Me,T(j,:),dof);
    
end

           
