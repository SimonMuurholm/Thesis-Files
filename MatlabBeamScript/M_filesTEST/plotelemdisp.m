function plotelemdisp(T,X,Ue,str,nonum)
%.........................................................
% plotelemdata:
%   Plots elements in topology matrix T with
%   coordinate matrix X. Uses linear line segment
%   between all data.
%
% Syntax:
%   plotelemdata(T,X,data,'linetype',nonum)
%
% Input:
%   T         :  element topology matrix.
%   X         :  node coordinate matrix.
%   Ue        :  displacements along the element
%   'linetype':  string defining linetype ('y-','c:',...).
%   nonum     :  if nonum = 1 -> node numbers will be plotted 
%
% Date:
%   Version 1.0    13.02.06
%.........................................................
% Node number offset
offset = 0.05*max(max(abs(X)));
%offset = 0.035*max(max(abs(X)));

% define line style and color
if nargin == 3
  str1 = 'yo';
  str2 = 'y-';
else
  if str(1) == ':' | str(1) == '-'  % check if line color is defined
    str1 = 'yo'; 
    str2 = ['y' str];
  else
    str1 = [str(1) 'o'];
    str2 = str;
  end
end

%%%%%
str3 = [str(1) ':'];
%%%%%


nnodes = size(T,2)-1;

% Plot 2D elements by calling function 'plot'
order = [1:nnodes];
% Plot nodes to scale geometry  
%
% plot(X(:,1),X(:,2),str1)
%
hold on

% Plot node numbers
if nonum==1
  for I=1:size(X,1)
    text(X(I,1)+offset,X(I,2)+offset,int2str(I))
  end
end

% Plot 2D elements
for j = 1:size(T,1)
%    
%  plot(X(T(j,order),1),X(T(j,order),2),str2)
   plot(X(T(j,order),1),X(T(j,order),2),str3)
%
end

% Plot data
ndata=size(Ue,2);
for j = 1:size(T,1)
  X1=X(T(j,1),:);
  a0=(X(T(j,2),:)-X(T(j,1),:));
  for k=1:ndata
    s=(k-1)/(ndata-1);
    Xd(k,:)=X1+a0*s+[Ue(j,k,1) Ue(j,k,2)];
  end
%%%  
%  plot(Xd(1:ndata,1),Xd(1:ndata,2),str2)
  plot(Xd(1:ndata,1),Xd(1:ndata,2),str2,'linewidth',1.0)
%%%
end

axis('equal')    % for use in ver. 4.0
axis('off') 

% enable the plot to be overwritten
hold off 
