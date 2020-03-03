function [x,v,a,t] = Newmark(K,C,M,x0,v0,dt,N,F,beta,gamma)
%.........................................................
%  NEWMARK
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
%    03.03.2015                 (MBNI)
%    Structural Engineering and Materials
%    Technical University of Denmark
%.........................................................

% ... Adjust default parameters 
n = length(x0);

% Use default parameters
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


% .. Modified mass matrix ...

MM = M + gamma*dt*C + beta*dt*dt*K; 
M1 = inv(MM);


% ... Initial values ...

t(1)   = 0.0;
x(:,1) =  x0;
v(:,1) =  v0;
a(:,1) =  M\(F(:,1) - C*v(:,1) - K*x(:,1));


% ... Time incrementation loop ...

for i = 1:N
   
  t(i+1) = t(i) + dt; 
  
  % Prediction step
  v(:,i+1) = v(:,i) + (1 - gamma)*dt*a(:,i);
  x(:,i+1) = x(:,i) + dt*v(:,i) + (0.5 - beta)*dt*dt*a(:,i);
  
  % Correction step
  a(:,i+1) = M1*(F(:,i+1) - C*v(:,i+1) - K*x(:,i+1));
  v(:,i+1) = v(:,i+1) + gamma*dt*a(:,i+1);
  x(:,i+1) = x(:,i+1) + beta*dt*dt*a(:,i+1);
  
end
