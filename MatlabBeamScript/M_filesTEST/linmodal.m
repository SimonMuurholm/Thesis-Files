function [x v t xs C mm km cm fm] = linmodal(K,C,M,x0,v0,omega0,S,t0,dt,N,f,alph,beta)
%.........................................................
% File: linmodal.m 
%    Function evaluating the forced dynamic response of 
%    a linear system with n degrees of freedom at equally 
%    spaced times using modal analysis.
%
%  INPUT
%    System matrices
%      K      : stiffness matrix.
%      M      : mass matrix.
%    Initial conditions
%      x0     : initial displacement.
%      v0     : initial velocity.
%    Vibration parameters
%      omega0 : angular frequency vector.
%      S      : Mode shape matrix.
%    Stepping parameters
%      t0     : initial time.
%      dt     : size of time step.
%      N      : number of time steps.
%    Load history
%      f      : vector of load amplitudes, optional.
%    Algorithm parameters
%      alph   : Parameter for stiffness proportional damping.
%      beta   : Parameter for mass proportional damping.
%
%  OUTPUT
%    Histories
%      x   : response history, x(n:N+1).
%      v   : velocity history, v(n:N+1).
%    Time
%      t   : discrete times, t(1:N+1).
%    Quasi-static correction
%      xs  : quasi-static response history.
%    System matrices
%      C   : damping matrix.
%    Modal matrices
%      mm   : modal mass matrix.
%      mm   : modal stiffness matrix.
%      cm   : modal damping matrix.
%      fm   : modal load history.
%
%
%  VERSION
%    01.02.2015
%    Structural Engineering and Materials
%    Technical University of Denmark
%.........................................................  

% mode shape matrix from frequency analysis
if ~exist('S')
    disp('... Run frequency analysis using linfreq.m first ...')
    return
end

% Dimension
n = size(S,2);

% Adjust load vector
if nargin < 11
    f = zeros(n,N+1);
else
    if size(f,2) < N+1
        f = [f, zeros(size(f,1),N+1-size(f,2))];
    end
end

%% ... Modal Analysis ...

% modal mass and stiffness
mm = S'*M*S;
km = S'*K*S;

% initial conditions
r0 = mm\S'*M*x0;
s0 = mm\S'*M*v0;

% modal load
fm = mm\(S'*f);

% damping matrix
if nargin == 13
C  = alph*K + beta*M;
end

cm = S'*C*S;
zeta = diag(cm)./( 2*diag(mm).*omega0 )

% modal time histories
r = zeros(n,N+1);
s = zeros(n,N+1);
for  j = 1:n
    
    % time response of mode j
    [rj,sj,t] = timestep2(omega0(j),zeta(j),r0(j),s0(j),t0,dt,N,fm(j,:));
    
    % store response in modal vector
    r(j,:) = rj;
    s(j,:) = sj;
    
end

% response history by superposition
x = S*r;
v = S*s;

% Quasi-static correction
xs = K\f;
for i = 1:n
    xs = xs - 1/km(i,i)*S(:,i)*S(:,i)'*f;
end


