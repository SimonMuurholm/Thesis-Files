function [xmax,xmin] = mmsort(x)
%--------------------------------------------------------------
%
% MMSort:
%   Identifies the maximum and minimum sequences, 'xmax' and 
%   'xmin' from a data record 'x'. Starts with a maximum and
%   truncates to equal number of maxima and minima.
%
% INPUT:    
%   x       :  Data sequence. 
%
% OUTPUT:
%   xmax    :  Sequence of maxima of x.
%   xmin    :  Sequence of minima of x.
%
% Version:
%   Department of Structural Engineering and Materials,
%   Technical University of Denmark.
%   12.06.1999.
%
%----------------------------------------------------------------

%----- Initialize variables etc. -------------------------------- 

N  = size(x,2);                  % Size of record

ymax = zeros(1,fix(N/2));        % Array for holding maxima
ymin = zeros(1,fix(N/2));        % Arrau for holding minima

%----- Identify 'internal' extremes -----------------------------

imax = 1;
imin = 1;

for i = 2:N-1
  if  (x(i-1) < x(i)) & (x(i) > x(i+1)) % local max
    ymax(imax) = x(i);
    imax = imax + 1; 
  elseif  (x(i-1) > x(i)) & (x(i) < x(i+1)) & (imax > 1) % local min
    ymin(imin) = x(i);
    imin = imin + 1; 
  end
end


%----- Removal of excess elements -------------------------------

M = min([imax,imin])-1;
xmax = ymax([1:M]);
xmin = ymin([1:M]);
