function [ic,iu] = constidx(C,dof,ndof)
%.........................................................
% constidx:
%   Evaluation of constrained/unconstrained index
%   sets
%
% Syntax:
%   [ic,iu] = constidx(C,dof)
%
% Input:
%   C    :  constraint definition array.
%   dof  :  degrees-of-freedom per node.
%   ndof :  total number of degrees-of-freedom.
%
% Output:
%   ic   :  index set of constrained dofs.
%   iu   :  index set of unconstrained dofs.
%
% Version 1.0    01.03.13
%.........................................................

% Constrained dofs
if size(C,1) >= 1  
    ic = (C(:,1)-1)*dof + C(:,2);
else
    ic    = zeros(1,size(C,1));
end

% Unconstrained dofs
iu = setdiff([1:ndof],ic);