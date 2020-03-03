function Uen = ubeam(T,X,H,U,Se,ns)
%.........................................................
% ubeam: 
%   Calculates displacements along the elements in a group 
%   of linear elastic beams.
%
% Syntax:
%   Ue = ubeam(T,X,G,u,n)
%
% Input:
%   T    :  topology matrix for bar elements.
%   X    :  node coordinate matrix. 
%   He   :  element properties [E A I (G As)]
%   U    :  element array of nodal displacements
%   Se   :  element section force array
%   ns   :  number of data points along the element
%
% Output:
%   Uen  :  displacements in elements
%
% Version 1.0    12.08.12
%.........................................................

% initialize element displacement array
Uen = zeros(size(T,1),ns,2);

% Loop over elements: calculate element displacements
for j = 1:size(T,1)
    
    % define element arrays
    Xe = X(T(j,1:2),:);
    He = H(T(j,3),:);
    
    % element displacement array
    Ue = U(T(j,1:2),:);
    
    % correction of nodal rotation by shear strain
    if size(He,2)>4 && He(4)>0
        Qe = Se(j,[1 end],2)'; % element shear force at nodes
        Ue(:,3) = Ue(:,3) - Qe/(He(4)*He(5));
    end
        
    % element displacement vector
    Ue = reshape(Ue',6,1);
    
    % Unit directional vector
    a0 = (Xe(2,:)-Xe(1,:))';
    L = sqrt(a0'*a0);
    n = a0/L;
    
    % Transformation metrix
    Ae = [n(1) n(2)  0    0    0   0
         -n(2) n(1)  0    0    0   0
           0    0    1    0    0   0
           0    0    0   n(1) n(2) 0
           0    0    0  -n(2) n(1) 0
           0    0    0    0    0   1];
    
    % Translational transformation matrix
    An = [n(1) n(2)
         -n(2) n(1)];
    
    % Displacements in local coordinates
    Ue = Ae*Ue;
    
    % Calculate global element displacements
    for i = 1:ns
        s = (i-1)/(ns-1);
        Ne = nebeam(Xe,s);
        Uen(j,i,:) = An'*Ne*Ue;
    end
end
