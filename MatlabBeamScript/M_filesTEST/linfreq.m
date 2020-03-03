function [omega0 S K M iu] = linfreq(X,T,H,C,dof,spring_dof_y,spring_dof_x,k_spring_y,k_spring_x)
%.........................................................
% File: linfreq.m 
%
%   Driver for vibration analysis of elastic frame
%   structure with beam elements.
%
%  INPUT
%      X      : Node coordinates.
%      T      : Topology.
%      H      : Material properties.
%      C      : Constraint matrix.
%    dof      : Nodal degrees of freedom.
%
%  OUTPUT
% omega0      : Angular frequency vector.
%      S      : Mode shape matrix (unnormalized).   
%      K      : System stiffness matrix (unconstrained).
%      M      : System mass matrix (unconstrained).
%      iu     : Vector of free dofs.
%
%  VERSION
%    01.02.2012
%    Structural Engineering and Materials
%    Technical University of Denmark
%.........................................................   

% Dimension
ndof = size(X,1)*dof;
[K,M] = deal(zeros(ndof));

% Identify constrained/unconstrained index sets
[ic,iu] = constidx(C,dof,ndof);

% Global stiffness matrix
K = kbeam(K,T,X,H,dof);

% Add spring
K(spring_dof_x,spring_dof_x) = K(spring_dof_x,spring_dof_x) + k_spring_x;
K(spring_dof_y,spring_dof_y) = K(spring_dof_y,spring_dof_y) + k_spring_y;

% Mass matrix
M = mbeam(M,T,X,H,dof);

% Vibration frequencies and modes ......
S = zeros(ndof,length(iu));
[S(iu,:),lam] = eig(K(iu,iu),M(iu,iu));

% Natural frequencies 
omega0 = sqrt(diag(lam));

% Arrange eigenvalus in ascending order
[omega0,idx] = sort(omega0);
S           = S(:,idx);


