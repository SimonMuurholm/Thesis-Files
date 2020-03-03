function [U Ue Se] = linstat(X,T,H,C,P,dof,spring_dof_y,spring_dof_x,k_spring_y,k_spring_x)
%.........................................................
% File: linstat.m 
%
%   Driver for linear analysis of elastic frame
%   structure with beam elements.
%
%  INPUT
%      X      : Node coordinates.
%      T      : Topology.
%      H      : Material properties.
%      C      : Prescribed displacements.
%      P      : Nodal loads.
%    dof      : Nodal degrees of freedom.
%   selem     : Elements to be used for super element generation.
%
% OUTPUT
%     Un      : Node displacement array.
%     Ue      : Element displacement array.
%     Se      : Forces in beam elements.
%
%  VERSION
%    01.02.2012
%    Structural Engineering and Materials
%    Technical University of Denmark
%.........................................................      
% Dimension
ndof = size(X,1)*dof;

% Identify constrained/unconstrained index sets
[ic,iu] = constidx(C,dof,ndof);

% Initialize matrices and vectors
u = zeros(ndof,1);
if size(C,2) == 3
    u(ic) = C(:,3);
end

if exist('selem')
    if ~exist('cdof')
        
        % Solve as superelement
        disp('Static Analysis by super element')
        Ks = selem.Ks;
        fs = selem.fs;
        Ts = selem.Ts;
        us = Ks\fs;
        u(iu) = Ts*us;
    
    elseif exist('cdof')
    
        disp('Static Analysis by Super Element Assembly')
        
        % Super element matrices
        Ks = selem.Ks;
        fs = selem.fs;
%         Ts = selem.Ts;
        
        % Nodal loads into load vector
        f = zeros(dof*size(X,1),1);
        f = loadnode(f,P,dof);
        
        % Global stiffness matrix
        K = kbeam(zeros(ndof),T,X,H,dof);

        % Assemble models
        KK = assemse(K,Ks,cdof);
        ff = assemse(f,fs,cdof);
        
        us = KK\ff;
        u = us;
%         u = [Ts*us(cdof); u];

    end
else

    disp('Static Analysis by Full Model')
    
    % Nodal loads into load vector
    f = zeros(dof*size(X,1),1);
    f = loadnode(f,P,dof);
    
    % Global stiffness matrix
    K = kbeam(zeros(ndof),T,X,H,dof);
    % Add spring 
    K(spring_dof_x,spring_dof_x) = K(spring_dof_x,spring_dof_x) + k_spring_x;
    K(spring_dof_y,spring_dof_y) = K(spring_dof_y,spring_dof_y) + k_spring_y;
    
    % Solve for unconstrained displacement components
    u(iu) = K(iu,iu)\(f(iu) - K(iu,ic)*u(ic));
end

% Nodal displacements
U = reshape(u,dof,size(u,1)/dof)';

% Element section forces
Se = sbeam(T,X,H,U,11);

% Element displacements
Ue = ubeam(T,X,H,U,Se,11);


