function K = kbeam(K,T,X,H,dof)
%.........................................................
% kbeam: 
%   Creates and assembles stiffness matrix
%   of elastic beam elements.
%
% Syntax:
%   K = kbeam(K,T,X,G)
%  
% Input:
%   K   :  existing global stiffness matrix.
%   T   :  topology matrix for beam elements.
%   X   :  nodal coordinate matrix
%   H   :  element property matrix H = [E A I (G As) me];
%   dof :  degrees of freedom per node
%
% Output:
%   K   :  new global stiffness matrix.
%
%.........................................................

% Global stiffness matrix by loop over elements
for j = 1:size(T,1)
    
    % Define element arrays
    Xe = X(T(j,1:2),:);       % element node coordinates
    He = H(T(j,3),:);         % element properties
    
    % element stiffness
    Ke = kebeam(Xe,He);
    
    % element stiffness into global format
    K = assem(K,Ke,T(j,:),dof);
    
end

           
