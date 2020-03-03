function [x,v,a,t,R] = GeneralizedAlpha(K,C,M,x0,v0,dt,N,F,beta,gamma,am,af)
%.........................................................
%  GENERALIZEDALPHA
%    Function evaluating the forced dynamic response of 
%    a linear system with n degrees of freedom at equally 
%    spaced times by direct collocation-based integration.
%
%  INPUT
%    System matrices
%      K      : stiffness matrix.
%      C      : damping matrix.
%      M      : mass matrix.
%    Initial conditions
%      x0     : initial displacement.
%      v0     : initial velocity.
%    Stepping parameters
%      dt     : size of time step.
%      N      : number of time steps.
%    Load history
%      F      : vector of load amplitudes, optional.
%    Algorithm parameters
%      beta   : Newmark parameter [0,0.5], default 1/4.
%      gamma  : Newmark parameter [0,1.0], default 1/2.
%      af     : Damping parameter of old inertia term.
%      am     : Damping paramerter on old force term.
%
%  OUTPUT
%    Histories
%      x   : response history, x(n:N+1).
%      v   : velocity history, v(n:N+1).
%      a   : acceleration history, a(n:N+1).
%    Time
%      t   : discrete times, t(1:N+1).
%
%  VERSION
%    05.09.2015
%    Structural Engineering and Materials
%    Technical University of Denmark
%.........................................................

% ... Adjust default parameters ...

n = length(x0);

% Use default parameters
if nargin < 12
  af = 0.0;
end 
if nargin < 11
  am = 0.0;
end 
if nargin < 10
  gamma = 0.5;
end 
if nargin < 9
  beta = 0.25; 
end   

% Adjust load vector
if nargin == 7
  F = zeros(1,N+1);
elseif nargin > 7 & size(F,2) < N+1
  F = [F, zeros(n,N+1-size(F,2))];
end


% ... Initialize output variables ...
x = zeros(n,N+1);
v = zeros(n,N+1);
a = zeros(n,N+1);
t = zeros(1,N+1);
R = zeros(3,N+1);

% .. Modified mass matrix ...
MM = (1-am)*M + (1-af)*(gamma*dt*C + beta*dt*dt*K); 
M1 = inv(MM);


% ... Initial values ...
t(1)   = 0.0;
x(:,1) =  x0;
v(:,1) =  v0;
a(:,1) =  M\(F(:,1) - C*v(:,1) - K*x(:,1));

% ... Time incrementation loop ...
for i = 1:N
   
  t(i+1) = t(i) + dt; 
  
  % Increment predictors
  dv = dt*a(:,i);
  dx = dt*v(:,i) + 1/2*dt^2*a(:,i);
  dF = F(:,i+1) - F(:,i);
  
  % Acceleration increment
  da = M1*( F(:,i) - ( M*a(:,i) + C*v(:,i) + K*x(:,i) ) ...
           + (1-af)*(dF - C*dv - K*dx)  );
  
  % State vector update
  a(:,i+1) = a(:,i) + da;
  v(:,i+1) = v(:,i) + dv + gamma*dt*da;
  x(:,i+1) = x(:,i) + dx + beta*dt^2*da;
  
  % Reactions
  RR = F(:,i+1) - ( M*a(:,i+1) + C*v(:,i+1));
  R(1,i+1) = sum(RR(1:3:end));
  R(2,i+1) = sum(RR(2:3:end));
  R(3,i+1) = sum(RR(3:3:end));
  
  if max(abs(da)) > 5*max(abs(da(:,1)))
      abs(da)
max(abs(da(:,1)))

      fprintf('Large increment causing instability: \n')
      fprintf('%1.2e,1.2e,1.2e,1.2e,1.2e,1.2e\n',diag(M(1:6,1:6)))
  end
end
