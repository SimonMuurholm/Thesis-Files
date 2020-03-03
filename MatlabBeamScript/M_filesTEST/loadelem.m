function f = loadelem(f,p,T,X,dof)
%***************************************************
% LoadElem:
%   Calculates and assembles equivalent element 
%   load vector in global coordinates.
%
% Input:
%   f   :  initial global load vector
%   p   :  element load array
%   T   :  topology matrix for beam elements
%   X   :  nodal coordinate matrix
%   dof :  degrees of freedom per node
%
% Output:
%   f   :  loads in global column vector
%
%***************************************************

% Equivalent nodal forces by loop over elements
for j = 1:size(p,1)
    
    % Element nodal coordinates
    Xe = X(T(p(j,1),1:2),:);
    
    % Global equivalent element load vector
    fe = febeam(p(j,2:3),Xe);
    
    % load vector into global format
    f = assem(f,fe,T(p(j,1),:),dof);
    
end