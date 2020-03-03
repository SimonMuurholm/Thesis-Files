%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Cantilever.m
%
%   Simple Cantilever for vibration analysis.
%
%   Data file to be read into memory before 
%   running the driver script 'DynaFrameV'.
%
% Version 1.0    01.09.15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear settings
close all
clear all
clc
 
%path(path,'..\M_files')
%...............................................
% Geometry
%...............................................
% Beam dimensions 
L = 80.0;
 
% Dimensions 
dof  = 3;       % Number of dof's per node
nno = 3;       % Number of nodes
nel = nno - 1;  % Number of elements
ndof = nno*dof; % Number of dof's
 
 
% Coordinates of nodes X = [x y],
X = [linspace(0,L,nno)' zeros(nno,1)];
 
% Topology matrix T = [node1 node2 propno],
T = [(1:(nno-1))'  (2:nno)'  ones(nel,1) ];
 
%...............................................
% Properties
%............................................... 
% Element property matrix H = [E  A  I (G As) m],
rho = 7.85e3;       % mass density
E   = 2.12e11;       % elastic modulus
nu = 0.3;           % Poisson ratio
G =  E/(2*(1+nu));  % shear modulus 
D_out = 6;          % Outer diameter
D_in = 5.8;         % Inner diameter
A_out = pi*(D_out/2)^2; % Outer area
A_in = pi*(D_in/2)^2;   % Inner area
A = A_out - A_in; % vertical member area
As = (A/2);      % vertical member shear area
I = (pi/64)*(D_out^4-D_in^4);      % vertical member moment of inertia
EI = E*I;           % Bending stiffness [Nm/m2]
EA = E*A;           % Axial stiffness [N]
EI_Orca = 1.71e9*1000;  % [N/m2]
EA_Orca = 393e6*1000;      % [N]
EI_relation = EI/EI_Orca;
EA_relation = EA/EA_Orca;
 
H = [ E  A  I G As rho*A];
  
 
%...............................................
% Static Loads
%............................................... 
% Concentrated loads P = [ node Px Py M ]
P = [ nno   0.000    1000e3    0.000 ]; % [N]
F = P(1,3);
F_Orca = F/1000; % [kN]
  
  
%...............................................
% Support conditions
%............................................... 
% Boundary conditions C = [node 'dof' u]  
C = [  1   1  0.0
       1   2  0.0
       1   3  0.0 ];
   
%...............................................
% Spring conditions
%...............................................   

spring_location = L; % Spring location along beam
dist    = abs(X(:,1) - spring_location); % Distance between nodes and spring
minDist = min(dist); % Smallest distance to node
node_num     = find(dist == minDist); % Finds the indices for the nearest node (row) in the X vector
spring_coordinate = X(node_num,1); % The coordinate of the node

angle = 0; % Cable angle [deg]
k_spring = 0; % [N/m] Spring stiffness
k_spring_x = k_spring*sin(angle*(pi/180)); % [N/m] Spring stiffness x direction
k_spring_y = k_spring*cos(angle*(pi/180)); % [N/m] Spring stiffness y direction
spring_node = node_num; % Node where spring should be applied
dof_local_x = 1; % Direction in which spring should work, 1=x, 2=y, 3=rotation
spring_dof_x = (spring_node - 1)*dof + dof_local_x; % Global dof for the applied spring
dof_local_y = 2; % Direction in which spring should work, 1=x, 2=y, 3=rotation
spring_dof_y = (spring_node - 1)*dof + dof_local_y; % Global dof for the applied spring   
 
 
%...............................................
% Parameters for dynamic analysis
%...............................................
% Initial conditions
x0 = zeros(ndof,1);
v0 = zeros(ndof,1);  
 
% time history
dt = 0.01; % [s]
N  = 3000; % [-]
t0 = 0;
Duration = dt*N; % [s]
 
% Rayleigh damping
alph = 0.00;
beta = 0.000;
 
 
 
%...............................................
% Plot parameters
%...............................................
% Axes used for geometry plots [Xmin Xmax Ymin Ymax]
PlotAxes = [-0.1*L 1*L  -0.25*L 0.25*L];
PlotFactorU = 1.1;
ploton = 0;


%...............................................
% Compare OrcaFlex response to Matlab response
% Load in OrcaFlex results
%...............................................

OpSystem = 'PC';
   switch(OpSystem)
   case 'Mac' 
    %% Import the data
[~, ~, raw] = xlsread('/Users/simonmuurholmhansen/Dropbox (WoodThilsted)/!E - Employees drive/SMH/Master Thesis/OrcaFlexValidation/ResultFiles/SineLoadResponseBeats.xlsx','Sheet1');
raw = raw(9:end,1:2);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

[~, ~, raw1] = xlsread('/Users/simonmuurholmhansen/Dropbox (WoodThilsted)/!E - Employees drive/SMH/Master Thesis/OrcaFlexValidation/ResultFiles/ResponseSTATIC.xlsx','Sheet1');
raw1 = raw1(9:end,1:2);
raw1(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw1)) = {''};

[~, ~, raw2] = xlsread('/Users/simonmuurholmhansen/Dropbox (WoodThilsted)/!E - Employees drive/SMH/Master Thesis/OrcaFlexValidation/ResultFiles/SpringTop.xlsx','Sheet1');
raw2 = raw2(9:end,1:2);
raw2(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw2)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

R1 = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw1); % Find non-numeric cells
raw1(R1) = {NaN}; % Replace non-numeric cells

R2 = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw2); % Find non-numeric cells
raw2(R2) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

data1 = reshape([raw1{:}],size(raw1));

data2 = reshape([raw2{:}],size(raw2));

%% Allocate imported array to column variable names
Time = data(:,1);
X1 = data(:,2);

TimeStatic = data1(:,1);
X1Static = data1(:,2);

TimeSpring = data2(:,1);
X1Spring = data2(:,2);

%% Clear temporary variables
clearvars data raw R;
clearvars data1 raw1 R1;
clearvars data2 raw2 R2;
 
   case 'PC' 
            %% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 2);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A19:B3019";

% Specify column names and types
opts.VariableNames = ["Time", "X1"];
opts.SelectedVariableNames = ["Time", "X1"];
opts.VariableTypes = ["double", "double"];

% Import the data

tb1 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\SineLoadResponseBeats_NEL_53.xlsx", opts, "UseExcel", false);
%tb1_400 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\SineLoadResponseBeats_NEL_400.xlsx", opts, "UseExcel", false);
tb2 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\ResponseSTATIC.xlsx", opts, "UseExcel", false);
tb3 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\SpringTop.xlsx", opts, "UseExcel", false);
tb4 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\SineLoadResponseLinGrowth.xlsx", opts, "UseExcel", false);
tb5 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\SineLoadResponse2freqSpring.xlsx", opts, "UseExcel", false);
tb6 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp10kN.xlsx", opts, "UseExcel", false);
tb7 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp100kN.xlsx", opts, "UseExcel", false);
tb8 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp1000kN.xlsx", opts, "UseExcel", false);
tb9 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp10000kN.xlsx", opts, "UseExcel", false);
tb10 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp100000kN.xlsx", opts, "UseExcel", false);
tb11 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp150000kN.xlsx", opts, "UseExcel", false);
tb12 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp200000kN.xlsx", opts, "UseExcel", false);
tb21 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp300000kN.xlsx", opts, "UseExcel", false);
tb13 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp400000kN.xlsx", opts, "UseExcel", false);
tb16 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp500000kN.xlsx", opts, "UseExcel", false);
tb17 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp600000kN.xlsx", opts, "UseExcel", false);
tb18 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp700000kN.xlsx", opts, "UseExcel", false);
tb19 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp800000kN.xlsx", opts, "UseExcel", false);
tb20 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp900000kN.xlsx", opts, "UseExcel", false);
tb14 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp1000000kN.xlsx", opts, "UseExcel", false);
tb22 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp1100000kN.xlsx", opts, "UseExcel", false);
tb23 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp1200000kN.xlsx", opts, "UseExcel", false);
tb24 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp1300000kN.xlsx", opts, "UseExcel", false);
tb25 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp1400000kN.xlsx", opts, "UseExcel", false);
tb15 = readtable("D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\OrcaFlexValidation\ResultFiles\LinCheck\SpringStatDisp1500000kN.xlsx", opts, "UseExcel", false);



%% Convert to output type

TimeStatic = tb2.Time;
X1Static = tb2.X1;

TimeSpring = tb3.Time;
X1Spring = tb3.X1;

% Time = tb4.Time;
% X1 = tb4.X1;

% This one is for the time dependent load
Time = tb5.Time;


%X1 = tb1.X1;
X1 = tb4.X1;


Time1 = tb6.Time;
X11 = tb6.X1;
X12 = tb7.X1;
X13 = tb8.X1;
X14 = tb9.X1;
X15 = tb10.X1;
X16 = tb11.X1;
X17 = tb12.X1;
X18 = tb13.X1;
X19 = tb14.X1;
X20 = tb15.X1;
X21 = tb16.X1;
X22 = tb17.X1;
X23 = tb18.X1;
X24 = tb19.X1;
X25 = tb20.X1;
X26 = tb21.X1;
X27 = tb22.X1;
X28 = tb23.X1;
X29 = tb24.X1;
X30 = tb25.X1;
%% Clear temporary variables
clear opts tbl tb2 tb3 tb4 tb5 tb6 tb7 tb8 tb9 tbl0 tbl1 tbl2 tbl3 tb14 tbl15 tbl16 tbl17 tbl18 tbl19 tbl20 tbl21 tbl22 tbl23 tbl24 tbl25
   end;
   
%umax_analytical = max(abs(x(fno*dof-1,:)));
umax_Orca = max(abs(X1-X1(1,1)));
 
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%...............................................
% Plot for check of geometric non linearity
%...............................................
loads = [0 10e3-k_spring_y*X11(11,1) 100e3-k_spring_y*X12(11,1) 1000e3-k_spring_y*X13(11,1) 10000e3-k_spring_y*X14(11,1) 100000e3-k_spring_y*X15(11,1) 200000e3-k_spring_y*X17(11,1) 300000e3-k_spring_y*X26(11,1) 400000e3-k_spring_y*X18(11,1) 500000e3-k_spring_y*X21(11,1) 600000e3-k_spring_y*X22(11,1) 700000e3-k_spring_y*X23(11,1) 800000e3-k_spring_y*X24(11,1) 900000e3-k_spring_y*X25(11,1) 1000000e3-k_spring_y*X19(11,1) 1100000e3-k_spring_y*X27(11,1) 1200000e3-k_spring_y*X28(11,1) 1300000e3-k_spring_y*X29(11,1) 1400000e3-k_spring_y*X30(11,1) 1500000e3-k_spring_y*X20(11,1)];
Displacements = [0 X11(11,1) X12(11,1) X13(11,1) X14(11,1) X15(11,1) X17(11,1) X26(11,1) X18(11,1) X21(11,1) X22(11,1) X23(11,1) X24(11,1) X25(11,1) X19(11,1) X27(11,1) X28(11,1) X29(11,1) X30(11,1) X20(11,1)];

%...............................................
% Linear Static
%...............................................

%AllPvals = [0 10e3 100e3 1000e3 10000e3 100000e3 1500000e3];
AllPvals = [0 1500000e3];
%AllPvals = loads;
datapoints = length(AllPvals);
savedUn = zeros(1,datapoints);
savedPvals = zeros(1,datapoints);
for i = 1:datapoints
    
    P(1,3) = AllPvals(1,i);
    
[Un Ue Se] = linstat(X,T,H,C,P,dof,spring_dof_y,spring_dof_x,k_spring_y,k_spring_x);
Un;
 
umax_analytical = P(1,3)*L^3/(3*E*I) + P(1,3)*L/(G*As) % Only valid for kspring = 0 (cantilever)
Un(end,2)

umax_analytical_new = (P(1,3)-k_spring_y*Un(end,2))*L^3/(3*E*I) % The 'resulting' force is used here
savedUn(1,i) = Un(end,2);
savedPvals(1,i) = P(1,3)-k_spring_y*Un(end,2);
 
end

%...............................................
% Frequency analysis
%...............................................
[omega0 S K M iu] = linfreq(X,T,H,C,dof,spring_dof_y,spring_dof_x,k_spring_y,k_spring_x);
C = alph*M + beta*K;
omega0(1);
Period1 = 2*pi/omega0(1);
%Period1_Orca = 1.05775; % First period
Period1_Orca = Period1; % First period
Freq1 = 1/Period1; % First frequency from Matlab script
Freq1_Orca = 1/Period1_Orca; % First frequency from OrcaFlex
omega1_analytical = 3.516*((E*I)/(rho*A*L^4))^0.5
 
%...............................................
% Dynamic Loads (Box or sine shaped load)
%...............................................
Load = 'Sine';
   switch(Load)
   case 'Box' 
        fd = zeros(ndof,N+1);
        fno = nno;
        Lper = 3; % Load period [s] 
        fd((fno*dof)-1,1:Lper/dt) = P(1,3)*ones(1,Lper/dt); % Place load on beam tip
                     
        if size(fd,2) < N+1
            fd = [fd, zeros(size(fd,1),N+1-length(fd))];
        end
 
   case 'Sine' 
      fd = zeros(ndof,N+1);
      fno = nno;
      %Lper = 3+dt; % Load period [s] 
      Lper = Duration+dt; % Load period [s]
      t_load = (0:dt:Lper-dt)';     % seconds/ time to be plotted
      %fc = Freq1*0.333;                     % Under freq, hertz / number of periods/s [1/T]
      fc = Freq1;                     % Linear growth, hertz / number of periods/s [1/T]
      fc_Orca = Freq1_Orca;
      %fc = Freq1*1.1;                     % Beats, hertz / number of periods/s [1/T]
      %fc_Orca = Freq1_Orca*1.1;                     % Beats, hertz / number of periods/s [1/T]
      %fc = Freq1*2;                     % Over freq, hertz / number of periods/s [1/T]
      %fc_Orca = Freq1_Orca*2;
      omega_c = fc*2*pi;                 % Angular frequency
      omega_t = sin(2*pi*fc*t_load); % omega = 2*pi/T so omega*t
      omega_t_Orca = sin(2*pi*fc_Orca*t_load); % omega = 2*pi/T so omega*t
      fd((fno*dof)-1,1:Lper/dt) = F*omega_t'.*ones(1,Lper/dt); % Place load on beam tip
      % Load placed on dof 2 of last node because it is y direction
      % ones(1,Lper/dt) inserts row of load values 
                     
      if size(fd,2) < N+1
      fd = [fd, zeros(size(fd,1),N+1-length(fd))];
   end;
   otherwise
      fprintf('Invalid load\n' );
   end
   
 
%...............................................
% Time integration
%...............................................
x0(iu,1) = K(iu,iu)\fd(iu,1);
lambda_inf = 0.80;
am = (2*lambda_inf-1)/(lambda_inf+1);
af = lambda_inf/(lambda_inf+1);
gamma = 1/2-am+af;
beta = 1/4*(1-am+af)^2;
 
% Time integration
[x(iu,:),v(iu,:),a(iu,:),t] = GeneralizedAlpha(K(iu,iu),C(iu,iu),M(iu,iu),x0(iu,1),v0(iu,1),dt,N,fd(iu,:),beta,gamma,am,af);
 
% 
% gamma = 0.5; beta = gamma/2;
% [x(iu,:),v(iu,:),a(iu,:),t] = Newmark(K(iu,iu),C(iu,iu),M(iu,iu),x0(iu,1),v0(iu,1),dt,N,fd(iu,:),beta,gamma);
 
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'DefaultAxesFontName','Times')
set(0,'DefaultAxesFontAngle','Normal')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontName','Times')
set(0,'DefaultTextFontAngle','Normal')
set(0,'DefaultTextFontSize',12)
 
 
%%...............................................
%% Plot for static analysis
%%...............................................
% Deformed geometry
figure(1); 
set(gcf,'windowstyle','docked')
subplot(221)
plotelem(T,X,'b:',1);
plotelemdisp(T,X,Ue*PlotFactorU,'b-',0);
axis(PlotAxes);
title('Initial and deformed geometry');
 
 
%...............................................
% Plot for frequency analysis
%...............................................
% Vibration mode
for mode = 1;
 
% Nodal vibration mode displacements
uv = S(:,mode)/max(abs(S(:,mode)));
Uvn = reshape(uv,dof,size(X,1))';
 
% Vibration mode section forces 
Sve = sbeam(T,X,H,Uvn,11);
 
% Vibration mode element displacements
Uve = ubeam(T,X,H,Uvn,Sve,11);
 
% Plot Vibration Modes
figure(2);
set(gcf,'windowstyle','docked')
subplot(2,2,mode)
plotelemdisp(T,X,Uve*0.0,'b:',0),hold on,
plotelemdisp(T,X,Uve,'b-',0)
title(['Mode ' num2str(mode), ...
      ', \omega_0 = ' num2str(omega0(mode),'%1.1f') ' rad/s'])
axis(PlotAxes);

figure(20);
set(gcf,'windowstyle','docked')
subplot(2,2,mode)
plotelemdisp(T,X,Uve*0.0,'b:',0),hold on,
plotelemdisp(T,X,Uve*3,'b-',0)
title(['Mode ' num2str(mode), ...
      ', \omega_0 = ' num2str(omega0(mode),'%1.1f') ' rad/s'])
axis(PlotAxes);
end
 

figure(21);
%set(gcf,'windowstyle','docked')
subplot(2,2,1)
plotelemdisp(T,X,Uve*0.0,'k:',0)
hold on
plotelemdisp(T,X,Uve*3,'k-',0)
%print('D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\DynamicAnalysisOfBuoyancyGuyedMonopiles\FIGURES\ValidationCases\ModeShape1MATLAB','-depsc')



%...............................................
% Plot for time domain analysis
%...............................................
figure(3);
set(gcf,'windowstyle','docked')
plot(t,x(fno*dof-1,:),'-'),hold on
%plot([0 t(end)],[Un(end,2) Un(end,2)],'r-')
xlabel('Time [s]')
ylabel('x_{tip} [m]')
 
%...............................................
% Plot load shape
%...............................................  
figure(4)
set(gcf,'windowstyle','docked')
plot(t,fd)
ylabel('Load [N]');
xlabel('Time [s]');
title('Load shape vs time');

%...............................................
% Extract sine load for OrcaFlex (when analysis run with sine load)
%...............................................  
Sine_load = (F*omega_t_Orca'.*ones(1,Lper/dt))/1000; % Size of load
Sine_Orca = zeros(length(Sine_load),2);
Sine_Orca(:,1) = t([1:length(Sine_Orca)]); % Time steps
Sine_Orca(:,2) = Sine_load; % Load steps

 

   
%...............................................
% Plot for time domain analysis response
%...............................................
% figure(5);
% %set(gcf,'windowstyle','docked')
% %subplot(221)
% p1 = plot(t,x(fno*dof-1,:),'-')
% hold on
% p2 = plot(Time,-X1,'-')
% %grid minor
% %p2 = plot(Time,X1_400,'-')
% %ylim([0 6])
% xlim([0 20])
% xlabel('t [s]','FontSize',16)
% ylabel('x(t) [m]','FontSize',16)
% ax = gca;
% ax.FontSize = 16;
% legend([p1 p2],'MATLAB','OrcaFlex','FontSize',16,'Location','Northeast')
% print('D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\DynamicAnalysisOfBuoyancyGuyedMonopiles\FIGURES\ValidationCases\ValOmegaBeats','-depsc')


figure(6);
%set(gcf,'windowstyle','docked')
%subplot(221)
p1 = plot(t,x(fno*dof-1,:),'-')
hold on
p2 = plot(Time,-X1,'-')
%grid minor
%ylim([0 6])
xlim([0 20])
xlabel('t [s]','FontSize',16)
ylabel('x(t) [m]','FontSize',16)
ax = gca;
ax.FontSize = 16;
legend([p1 p2],'MATLAB','OrcaFlex','FontSize',16,'Location','Northwest')
print('D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\DynamicAnalysisOfBuoyancyGuyedMonopiles\FIGURES\ValidationCases\ValOmegaLinGrow','-depsc')




%...............................................
% Plot for check of geometric non linearity
%...............................................
%loads = [0 10e3-k_spring_y*X11(11,1) 100e3-k_spring_y*X12(11,1) 1000e3-k_spring_y*X13(11,1) 10000e3-k_spring_y*X14(11,1) 100000e3-k_spring_y*X15(11,1) 150000e3-k_spring_y*X16(11,1) 200000e3-k_spring_y*X17(11,1) 400000e3-k_spring_y*X18(11,1) 1000000e3-k_spring_y*X19(11,1) 1500000e3-k_spring_y*X20(11,1)];
%loads = [0 10 100 1000 10000 100000 150000 200000 400000];

% figure(9);
% set(gcf,'windowstyle','docked')
% subplot(221)
% p1 = plot(loads/1000,Displacements,'-o'),hold on
% p2 = plot(savedPvals/1000,savedUn,'-o')
% %p1 = plot(Displacements,loads/1000,'-o'),hold on
% %p2 = plot(savedUn,savedPvals/1000,'-o')
% %plot([0 loads(end)/1000],[0 Displacements(end)],'r-') % Line between first and last point
% %p2 = plot(Time,X1,'-')
% %ylim([0 6])
% %xlim([0 1000])
% %ylim([0 20])
% xlabel('Resulting load [kN]','FontSize',12)
% ylabel('x_{tip} (x-direction) [m]','FontSize',12)
% ax = gca;
% ax.FontSize = 12;
% legend([p1 p2],'OrcaFlex','Matlab','Location','Northwest')
% %print('/Users/simonmuurholmhansen/Dropbox (WoodThilsted)/!E - Employees drive/SMH/Master Thesis/DynamicAnalysisOfBuoyancyGuyedMonopiles/FIGURES/ValidationCases/Test','-depsc')
% print('D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\DynamicAnalysisOfBuoyancyGuyedMonopiles\FIGURES\ValidationCases\NonLinCheck','-depsc')


% figure(10);
% set(gcf,'windowstyle','docked')
% subplot(221)
% p1 = plot(loads/1000,Displacements,'-o'),hold on
% p2 = plot(savedPvals/1000,savedUn,'-o')
% %p1 = plot(Displacements,loads/1000,'-o'),hold on
% %p2 = plot(savedUn,savedPvals/1000,'-o')
% %plot([0 loads(end)/1000],[0 Displacements(end)],'r-') % Line between first and last point
% %p2 = plot(Time,X1,'-')
% %ylim([0 6])
% xlim([0 1000])
% %ylim([0 20])
% xlabel('Resulting load [kN]','FontSize',12)
% ylabel('x_{tip} (x-direction) [m]','FontSize',12)
% ax = gca;
% ax.FontSize = 12;
% legend([p1 p2],'OrcaFlex','Matlab','Location','Northwest')
% %print('/Users/simonmuurholmhansen/Dropbox (WoodThilsted)/!E - Employees drive/SMH/Master Thesis/DynamicAnalysisOfBuoyancyGuyedMonopiles/FIGURES/ValidationCases/Test','-depsc')
% print('D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\DynamicAnalysisOfBuoyancyGuyedMonopiles\FIGURES\ValidationCases\NonLinCheckZoom','-depsc')



% figure(11);
% set(gcf,'windowstyle','docked')
% subplot(221)
% %p1 = plot(loads/1000,Displacements,'-o'),hold on
% %p2 = plot(savedPvals/1000,savedUn,'-o')
% p1 = plot(Displacements/L,(loads*L^2)/(EI_Orca),'-x'),hold on
% p2 = plot(savedUn/L,(savedPvals*L^2)/(EI),'-x')
% %plot([0 loads(end)/1000],[0 Displacements(end)],'r-') % Line between first and last point
% %p2 = plot(Time,X1,'-')
% %ylim([0 6])
% %xlim([0 20])
% %ylim([0 20])
% xlabel('u_1/L','FontSize',12)
% ylabel('PL^2/EI','FontSize',12)
% ax = gca;
% ax.FontSize = 12;
% legend([p1 p2],'OrcaFlex','Matlab','Location','Northwest')
% %print('/Users/simonmuurholmhansen/Dropbox (WoodThilsted)/!E - Employees drive/SMH/Master Thesis/DynamicAnalysisOfBuoyancyGuyedMonopiles/FIGURES/ValidationCases/NonLinCheck','-depsc')
% print('D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\DynamicAnalysisOfBuoyancyGuyedMonopiles\FIGURES\ValidationCases\NonLinCheckNorm','-depsc')


%...............................................
% Calculate dynamic amplification Matlab
%...............................................

N_points = length(-x(fno*dof-1,:));
omega_c_vec = [0*omega_c:9*omega_c/(N_points-1):9*omega_c];
r_vec = omega_c_vec./omega0(1);
DynAmpl_vec = 1./(sqrt((1-r_vec.^2).^2));
% ADD HERE THE DAF FOR DAMPED RESPONSE

t_norm = t./Period1;

%...............................................
% Normalize response Matlab
%...............................................

F_vec = zeros(ndof,1);
F_vec((fno*dof)-1,1) = F;

S(:,1) = S(:,1)/max(abs(S(:,1))); % Normalize first modeshape vector
f_val = S(:,1)'*F_vec;
k_val = S(:,1)'*K*S(:,1);
f_div_k = f_val/k_val;
k_div_f = k_val/f_val;
Response = -x(fno*dof-1,:);
Norm_response = -Response*k_div_f;

%...............................................
% Normalize response OrcaFlex
%...............................................

Period1Orca = 1.058; %[s] First period
t_normOrca = Time./Period1Orca; % Normalised time [-]
Eqv_Stat = X1Static(1)*(-F); %[m] Equivalent static displacement from unit load 1 [N] multiplied with load F 
Norm_responseOrca = (X1*(-1))/Eqv_Stat; % Dynamic response normalised with equivalent static load [-]


% figure(6);
% set(gcf,'windowstyle','docked')
% %subplot(2,2,[3,4]);
% subplot(2,2,1);
% p1 = plot(t_norm,Norm_response,'-')
% hold on
% p2 = plot(t_normOrca,Norm_responseOrca,'-') 
% %ylim([0 6])
% xlim([0 20])
% xlabel('t/T_1 [-]','FontSize',12)
% ylabel('x(t) k/F [-]','FontSize',12)
% ax = gca;
% ax.FontSize = 12;
% legend([p1 p2],'Matlab','OrcaFlex')
% print('/Users/simonmuurholmhansen/Dropbox (WoodThilsted)/!E - Employees drive/SMH/Master Thesis/DynamicAnalysisOfBuoyancyGuyedMonopiles/FIGURES/ValidationCases/ValCase1Normalized','-depsc')


% figure(6);
% set(gcf,'windowstyle','docked')
% %subplot(2,2,[3,4]);
% %subplot(2,2,1);
% p1 = plot(t_norm,Norm_response,'-')
% hold on
% p2 = plot(t_normOrca,Norm_responseOrca,'-') 
% %ylim([0 6])
% xlim([0 20])
% xlabel('t/T_1 [-]','FontSize',12)
% ylabel('x(t) k/F [-]','FontSize',12)
% ax = gca;
% ax.FontSize = 12;
% legend([p1 p2],'Matlab','OrcaFlex')
%print('/Users/simonmuurholmhansen/Dropbox (WoodThilsted)/!E - Employees drive/SMH/Master Thesis/DynamicAnalysisOfBuoyancyGuyedMonopiles/FIGURES/ValidationCases/ValCase1NormalizedTEST','-depsc')
%print('D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\DynamicAnalysisOfBuoyancyGuyedMonopiles\FIGURES\ValidationCases\ValCase1NormalizedLinGrowth','-depsc')


%...............................................
% Plot for dynamic amplification 
%...............................................

figure(7);
%set(gcf,'windowstyle','docked')
%subplot(221)
plot(r_vec,DynAmpl_vec,'k-')
grid minor
ylim([0 6])
xlim([0 3])
xlabel('f_{ex}/f_0 [-]','FontSize',16)
ylabel('DAF [-]','FontSize',16)
ax = gca;
ax.FontSize = 16;
%legend([p1 p2],'MATLAB','OrcaFlex')
print('D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\DynamicAnalysisOfBuoyancyGuyedMonopiles\FIGURES\ValidationCases\DynamicAmplification','-depsc')

%...............................................
% Plot for spring comparisson Matlab and OrcaFlex
%...............................................
% 
% figure(8);
% set(gcf,'windowstyle','docked')
% %subplot(221)
% p1 = plot(t,-x(fno*dof-1,:),'-'),hold on
% %plot([0 t(end)],[Un(end,2) Un(end,2)],'r-')
% p2 = plot(TimeSpring,X1Spring*(-1),'-')
% %ylim([0 6])
% xlim([0 20])
% xlabel('Time [s]')
% ylabel('x_{tip} [m]')
% legend([p1 p2],'Matlab','OrcaFlex')


