function A = assemse(A,Ae,cdof)
%.........................................................
% Assem: 
%   Assembles system matrix or vector by adding 
%   contributions from super elements to an 
%   existing global matrix. The system matrix may be 
%   square like the stiffness mass or conductivity 
%   matrix, or may be a one-dimensional vector like 
%   the load vector.
%
% Syntax:
%   A = assmk(A,Ae,cdof)
%
% Input:
%   A  :  global matrix.
%   Ae :  element matrix.
%  cdof:  connection dofs (global).
%
% Output:
%   A  :  updated global matrix. 
%
% Date:
%   Version 1.0    03.09.15
%.........................................................

% Number of dofs
ndof = size(A,1);

% Number of dof in Super Element
nsedof = size(Ae,1);

% Define global address vector for element dofs
ig = [cdof, ndof+1:ndof+nsedof-length(cdof)];

% Add element matrix/vector to global matrix/vector
if size(A,2) ~= size(A,1)
    AA = zeros(ndof+nsedof-length(cdof),1);
    AA(1:ndof,1) = A;
    A(ig) = AA(ig) + Ae;
else
    AA = zeros(ndof+nsedof-length(cdof));
    AA(1:ndof,1:ndof) = A;
    A(ig,ig) = AA(ig,ig) + Ae;
end
     
