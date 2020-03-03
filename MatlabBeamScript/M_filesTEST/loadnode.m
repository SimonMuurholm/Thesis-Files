function f = loadnode(f,P,dof)
%.........................................................
% LoadNode: 
%   Sets the nodal loads in global vector form. 
%
% Input:
%   f   :  initial global load vector
%   P   :  nodal load array
%   dof :  degrees of freedom per node
%
% Output:
%   f   :  loads in global column vector
%
% Date:
%   Version 1.0    04.05.95
%   Version 1.1    17.02.05
%.........................................................

% Put loads into global load vector.
for i = 1:size(P,1)
    j = (P(i,1)-1)*dof;
    f(j+1:j+dof) = P(i,1+1:1+dof)';
end