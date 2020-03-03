function A = assem(A,Ae,Te,dof)
%.........................................................
% Assem: 
%   Assembles system matrix or vector by adding 
%   element contributions from elements to an 
%   existing global matrix. The system matrix may be 
%   square like the stiffness mass or conductivity 
%   matrix, or may be a one-dimensional vector like 
%   the load vector.
%
% Syntax:
%   A = assmk(A,Ae,Te)
%   A = assmk(A,Ae,Te,dof)
%
% Input:
%   A  :  global matrix.
%   Ae :  element matrix.
%   Te :  element topology vector.
%   dof:  degrees of freedom per node.
%
% Output:
%   A  :  updated global matrix. 
%
% Date:
%   Version 1.0    10.08.12
%.........................................................

% Number of element nodes
enodes = size(Te,2) - 1; 

% Define global address vector for element dofs
ig = zeros(1,enodes*dof);
for i = 1:enodes
  ig([(i-1)*dof+1:i*dof]) = [(Te(i)-1)*dof+1:Te(i)*dof]; 
end

% Add element matrix/vector to global matrix/vector
if size(A,2) ~= size(A,1) 
  A(ig) = A(ig) + Ae;
else
  A(ig,ig) = A(ig,ig) + Ae;
end
     
