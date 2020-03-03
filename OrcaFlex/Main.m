% This is the main file which runs the provided script 
% and produces results
clc
clear all
close all

addpath('C:\Program Files (x86)\Orcina\OrcaFlex\11.0\OrcFxAPI\MATLAB')

filename = 'D:\Dropbox (WoodThilsted)\!E - Employees drive\SMH\Master Thesis\ScriptedOrcaFlex\AllModels\OutputFiles\OrcaFlexModel';
tic
[NumLines,MSL,LoadAngle,lvals,Lvals,Time,umaxMP,l,L,EndForceResults,uMP,... 
    EndForce,umaxMPAllGeometry,EndForceMaxAllGeometryChain1,EndForceMaxAllGeometryMP,... 
    MomentMPAllGeometry,VonMisesMPAllGeometry,ShearMPAllGeometry,umaxTowerAllGeometry,...
    EndForceMaxAllGeometryChain4,etaAllGeometry,EffectiveTensionMPAllGeometry,...
    WallTensionMPAllGeometry,ExternalPressureMPAllGeometry,TotalFoundationLength,...
    EndForceMaxAllGeometryChain4EndA,EndForceTimeHistory,VonMisesTimeAllGeomtery,...
    uTowerAllGeometry] = BuildModel(filename);

[GeometryPlot,ch,msl,width,height,yint] = DrawGeometry(); % Draws geometry 
fprintf('Analysis finished \n')
toc


%%%%% Find frequency of waves %%%%%%%%%

eta_ts(:,1) = Time(1,:);
eta_ts(:,2) = etaAllGeometry(1,:);

% Finds the number of data-set in the excel-file
points = size(eta_ts(:,1));

% We go through the time series to find the downcrossings
% Therefore we set-up a for-loop 
kk=0;
for ii=1+1:points
    % Test if our water level no. ii is negative, while the previous (ii-1) was possitive.
    % If the if-statement below is true a down-crossing point is found
    if (eta_ts(ii,2)<0 && eta_ts(ii-1,2) >0)  
       kk=kk+1;        % kk is the number of down-crossings. As we have found another down-crossing we add 1 to kk
       zero_c_n(kk) = ii; % variable used to remember at which data-set down-crossings are found.
       % linear interpolation finds down-crossing in next three lines
       w1 = (0-eta_ts(ii,2))/(eta_ts(ii-1,2)-eta_ts(ii,2));
       w2 = 1-w1;
       time(kk)     = w2*(eta_ts(ii,1)-eta_ts(ii-1,1)) +eta_ts(ii-1,1);
       % The level is set equal to 0. 
       zc_z(kk) = 0.0;
    end
end
%  Set kmax = kk, to save the maximum number of down-crossings
kmax = kk;

% Plot the time-series with downcrossings:
figure(100)
    plot( eta_ts(:,1),eta_ts(:,2),'k','LineWidth',1)
    hold on  
    plot( time,zc_z,'or','LineWidth',1,'MarkerSize',4)
    grid minor
    ylim([-7 6])
    xlabel('Time [s]','FontSize',16)
    ylabel('Sea elevation [m]','FontSize',16)
    ax = gca;
    ax.FontSize = 16;
    str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\zerodown',NumLines);
    print(str,'-dpng')
    
figure(107)
    plot( eta_ts(:,1),eta_ts(:,2),'k','LineWidth',1)
    hold on
    plot( time,zc_z,'or','LineWidth',1,'MarkerSize',4)
    grid minor
    xlim([190 340])
    ylim([-7 6])
    xlabel('Time [s]','FontSize',16)
    ylabel('Sea elevation [m]','FontSize',16)
    ax = gca;
    ax.FontSize = 16;
    str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\zerodown2',NumLines);
    print(str,'-dpng')
% Time series analyses in the time domain:

    f1_25m = 0.158;
    f2_25m = 0.759;
    zeta = 0.008;
% Find wave heights and wave periods
for i=1:length(time)-1
    WavePeriod(i) = time(i+1)-time(i);
    WaveHeight(i) = max(eta_ts(zero_c_n(i):zero_c_n(i+1),2))+abs(min(eta_ts(zero_c_n(i):zero_c_n(i+1),2)));
    WaveStat(i,1)=WavePeriod(i);
    WaveStat(i,2)=WaveHeight(i); 
    WaveStat(i,3)=i; %Wave number
    
    fw(i) = 1/WavePeriod(i);   % Wave frequencies
    fw_save(i,1) = 1/WavePeriod(i);   % Wave frequencies
    r1(i) = fw(i)/f1_25m;
    r2(i) = fw(i)/f2_25m;
    
    DAF1(i) = ( 1/ ( sqrt((1-(r1(i))^2)^2) + (2*zeta*r1(i))^2));
    DAF2(i) = ( 1/ ( sqrt((1-(r2(i))^2)^2) + (2*zeta*r2(i))^2));
end


figure(101)
plot(time(1:end-1),1./WaveStat(:,1),'ok-')
grid minor
xlabel('Time [s]','FontSize',16)
ylabel('Wave frequency [Hz]','FontSize',16)
ax = gca;
ax.FontSize = 16;
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\WaveFrequencies',NumLines);
print(str,'-dpng')

figure(102)
plot(time(1:end-1),1./WaveStat(:,1),'ok-')
grid minor
xlabel('Time [s]','FontSize',16)
ylabel('Wave frequency [Hz]','FontSize',16)
xlim([190 340])
ax = gca;
ax.FontSize = 16;
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\WaveFrequenciesZoom',NumLines);
print(str,'-dpng')



figure(103)
plot(time(1:end-1),DAF1,'or-')
hold on 
plot(time(1:end-1),DAF2,'ob-')
grid minor
xlabel('Time [s]','FontSize',16)
ylabel('DAF [-]','FontSize',16)
ax = gca;
ax.FontSize = 16;
legend({'Mode 1','Mode 2'},'FontSize',16,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\WaveDAF',NumLines);
print(str,'-dpng')

figure(104)
plot(time(1:end-1),DAF1,'or-')
hold on 
plot(time(1:end-1),DAF2,'ob-')
grid minor
xlabel('Time [s]','FontSize',16)
ylabel('DAF [-]','FontSize',16)
xlim([190 340])
ax = gca;
ax.FontSize = 16;
legend({'Mode 1','Mode 2'},'FontSize',16,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\WaveDAF2',NumLines);
print(str,'-dpng')

figure(108)
plot(time(1:end-1),DAF1,'ok-')
hold on 
plot(time(1:end-1),DAF2,'ok--')
grid minor
xlabel('Time [s]','FontSize',16)
ylabel('DAF [-]','FontSize',16)
xlim([190 340])
ax = gca;
ax.FontSize = 16;
legend({'Mode 1','Mode 2'},'FontSize',16,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\WaveDAF2v2',NumLines);
print(str,'-dpng')


figure(105)
plot(time(1:end-1),r1,'or-')
hold on 
plot(time(1:end-1),r2,'ob-')
grid minor
xlabel('Time [s]','FontSize',16)
ylabel('f_w/f_n [-]','FontSize',16)
ax = gca;
ax.FontSize = 16;
legend({'Mode 1','Mode 2'},'FontSize',16,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\WaveRrelation',NumLines);
print(str,'-dpng')


figure(106)
plot(time(1:end-1),r1,'or-')
hold on 
plot(time(1:end-1),r2,'ob-')
grid minor
xlabel('Time [s]','FontSize',16)
ylabel('f_w/f_n [-]','FontSize',16)
xlim([190 340])
ylim([0 1.3])
ax = gca;
ax.FontSize = 16;
legend({'Mode 1','Mode 2'},'FontSize',16,'Location','Northeast')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\WaveRrelation2',NumLines);
print(str,'-dpng')

figure(110)
plot(time(1:end-1),r1,'ok-')
hold on 
plot(time(1:end-1),r2,'ok--')
grid minor
xlabel('Time [s]','FontSize',16)
ylabel('f_w/f_n [-]','FontSize',16)
xlim([190 340])
ylim([0 1.3])
ax = gca;
ax.FontSize = 16;
legend({'Mode 1','Mode 2'},'FontSize',16,'Location','Northeast')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\WaveRrelation2v2',NumLines);
print(str,'-dpng')



%%%%%%%%%%%%% OrcaFlex Results %%%%%%%%%%%%%%

x0=2;
y0=2;
widthNEW=24;
heightNEW=9.5;

figure(1)
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
y1 = [umaxMPAllGeometry(1,:)*1000];
x1 = [25];
y2 = [umaxMPAllGeometry(2,:)*1000];
x2 = [30];
y3 = [umaxMPAllGeometry(3,:)*1000];
x3 = [40];
y4 = [umaxMPAllGeometry(4,:)*1000];
x4 = [50];
y5 = [umaxMPAllGeometry(5,:)*1000];
x5 = [60];
y6 = [umaxMPAllGeometry(6,:)*1000];
x6 = [70];
y7 = [umaxMPAllGeometry(7,:)*1000];
x7 = [80];
y8 = [umaxMPAllGeometry(8,:)*1000];
x8 = [90];
bar(x1,y1,3.5)
hold on
bar(x2,y2,3.5)
bar(x3,y3,3.5)
bar(x4,y4,3.5)
bar(x5,y5,3.5)
bar(x6,y6,3.5)
bar(x7,y7,3.5)
bar(x8,y8,3.5)
bar(113,0)
grid minor
xlabel('l [m]','FontSize',15)
ylabel('Max displacement [mm]','FontSize',15)
ax = gca;
ax.FontSize = 15;
xticks([25 30 40 50 60 70 80 90 150])
xticklabels({'25','30','40','50','60','70','80','90',''})
legend({'L/l = 1.985','L/l = 1.745','L/l = 1.465','L/l = 1.316','L/l = 1.227','L/l = 1.170','L/l = 1.132','L/l = 1.105'},'FontSize',15,'Location','Northeast')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\umax_MP_HingedBot_TEST1',NumLines);
print(str,'-dpng')


figure(2)
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
y1 = [EndForceMaxAllGeometryChain4(1,:)/1000];
x1 = [25];
y2 = [EndForceMaxAllGeometryChain4(2,:)/1000];
x2 = [30];
y3 = [EndForceMaxAllGeometryChain4(3,:)/1000];
x3 = [40];
y4 = [EndForceMaxAllGeometryChain4(4,:)/1000];
x4 = [50];
y5 = [EndForceMaxAllGeometryChain4(5,:)/1000];
x5 = [60];
y6 = [EndForceMaxAllGeometryChain4(6,:)/1000];
x6 = [70];
y7 = [EndForceMaxAllGeometryChain4(7,:)/1000];
x7 = [80];
y8 = [EndForceMaxAllGeometryChain4(8,:)/1000];
x8 = [90];
bar(x1,y1,3.5)
hold on
bar(x2,y2,3.5)
bar(x3,y3,3.5)
bar(x4,y4,3.5)
bar(x5,y5,3.5)
bar(x6,y6,3.5)
bar(x7,y7,3.5)
bar(x8,y8,3.5)
bar(113,0)
yl = yline(19.2,'--k','Load limit','LineWidth',3,'FontSize',15);
yl.LabelHorizontalAlignment = 'center';
grid minor
ylim([0 25])
xlabel('l [m]','FontSize',15)
ylabel('Max chain force [MN]','FontSize',15)
ax = gca;
ax.FontSize = 15;
xticks([25 30 40 50 60 70 80 90 120])
xticklabels({'25','30','40','50','60','70','80','90',''})
legend({'L/l = 1.985','L/l = 1.745','L/l = 1.465','L/l = 1.316','L/l = 1.227','L/l = 1.170','L/l = 1.132','L/l = 1.105'},'FontSize',15,'Location','Northeast')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\Chain4_MaxEndForce_EndB_HingedBot_TEST1',NumLines);
print(str,'-dpng')



figure(3)
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
y1 = [-EndForceMaxAllGeometryMP(1,:)/1000];
x1 = [25];
y2 = [-EndForceMaxAllGeometryMP(2,:)/1000];
x2 = [30];
y3 = [-EndForceMaxAllGeometryMP(3,:)/1000];
x3 = [40];
y4 = [-EndForceMaxAllGeometryMP(4,:)/1000];
x4 = [50];
y5 = [-EndForceMaxAllGeometryMP(5,:)/1000];
x5 = [60];
y6 = [-EndForceMaxAllGeometryMP(6,:)/1000];
x6 = [70];
y7 = [-EndForceMaxAllGeometryMP(7,:)/1000];
x7 = [80];
y8 = [-EndForceMaxAllGeometryMP(8,:)/1000];
x8 = [90];
bar(x1,y1,3.5)
hold on
bar(x2,y2,3.5)
bar(x3,y3,3.5)
bar(x4,y4,3.5)
bar(x5,y5,3.5)
bar(x6,y6,3.5)
bar(x7,y7,3.5)
bar(x8,y8,3.5)
bar(113,0)
grid minor
xlabel('l [m]','FontSize',15)
ylabel('Max reaction force [MN]','FontSize',15)
ax = gca;
ax.FontSize = 15;
xticks([25 30 40 50 60 70 80 90 120])
xticklabels({'25','30','40','50','60','70','80','90',''})
legend({'L/l = 1.985','L/l = 1.745','L/l = 1.465','L/l = 1.316','L/l = 1.227','L/l = 1.170','L/l = 1.132','L/l = 1.105'},'FontSize',15,'Location','Northeast')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_MaxEndForce_EndA_HingedBot_TEST1',NumLines);
print(str,'-dpng')

figure(7)
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
y1 = [umaxTowerAllGeometry(1,:)*1000];
x1 = [25];
y2 = [umaxTowerAllGeometry(2,:)*1000];
x2 = [30];
y3 = [umaxTowerAllGeometry(3,:)*1000];
x3 = [40];
y4 = [umaxTowerAllGeometry(4,:)*1000];
x4 = [50];
y5 = [umaxTowerAllGeometry(5,:)*1000];
x5 = [60];
y6 = [umaxTowerAllGeometry(6,:)*1000];
x6 = [70];
y7 = [umaxTowerAllGeometry(7,:)*1000];
x7 = [80];
y8 = [umaxTowerAllGeometry(8,:)*1000];
x8 = [90];
bar(x1,y1,3.5)
hold on
bar(x2,y2,3.5)
bar(x3,y3,3.5)
bar(x4,y4,3.5)
bar(x5,y5,3.5)
bar(x6,y6,3.5)
bar(x7,y7,3.5)
bar(x8,y8,3.5)
bar(113,0)
grid minor
xlabel('l [m]','FontSize',15)
ylabel('Displacement [mm]','FontSize',15)
ax = gca;
ax.FontSize = 15;
xticks([25 30 40 50 60 70 80 90 120])
xticklabels({'25','30','40','50','60','70','80','90',''})
legend({'L/l = 1.985','L/l = 1.745','L/l = 1.465','L/l = 1.316','L/l = 1.227','L/l = 1.170','L/l = 1.132','L/l = 1.105'},'FontSize',15,'Location','Northeast')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\umax_Tower_HingedBot_TEST1',NumLines);
print(str,'-dpng')


figure(8);
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
plot(Time,etaAllGeometry(1,:),'-k');
hold on 
grid minor
ylim([-7 6])
xlabel('Time [s]','FontSize',15)
ylabel('Sea elevation [m]','FontSize',15)
ax = gca;
ax.FontSize = 15;
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\WaveElevation_TEST1',NumLines);
print(str,'-dpng')

figure(30);
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
p1 = plot(Time,EndForceTimeHistory(1,:)/1000,'-');
hold on 
p2 = plot(Time,EndForceTimeHistory(2,:)/1000,'-');
p3 = plot(Time,EndForceTimeHistory(3,:)/1000,'-');
p4 = plot(Time,EndForceTimeHistory(4,:)/1000,'-');
p5 = plot(Time,EndForceTimeHistory(5,:)/1000,'-');
p6 = plot(Time,EndForceTimeHistory(6,:)/1000,'-');
p7 = plot(Time,EndForceTimeHistory(7,:)/1000,'-');
p8 = plot(Time,EndForceTimeHistory(8,:)/1000,'-');
grid minor
xlabel('Time [s]','FontSize',15)
ylabel('Fairlead force [MN]','FontSize',15)
xlim([0 max(Time*1.2)]);
ax = gca;
ax.FontSize = 15;
legend({'L/l = 1.985','L/l = 1.745','L/l = 1.465','L/l = 1.316','L/l = 1.227','L/l = 1.170','L/l = 1.132','L/l = 1.105'},'FontSize',15,'Location','Northeast')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\EndForceChain4TimeHistory_TEST1',NumLines);
print(str,'-dpng')



figure(33);
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
p1 = plot(Time,EndForceTimeHistory(4,:)/1000,'-');
hold on 
p2 = plot(Time,EndForceTimeHistory(1,:)/1000,'-');
grid minor
xlabel('Time [s]','FontSize',15)
ylabel('Fairlead force [MN]','FontSize',15)
xlim([0 max(Time)]);
ax = gca;
ax.FontSize = 15;
legend({'l = 50 m','l = 25 m'},'FontSize',15,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\EndForceChain4TimeHistorySingleChainLowHigh_TEST1',NumLines);
print(str,'-dpng')

figure(56);
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
p1 = plot(Time,EndForceTimeHistory(4,:)/1000,'--k');
hold on 
p2 = plot(Time,EndForceTimeHistory(1,:)/1000,'-k');
grid minor
xlabel('Time [s]','FontSize',15)
ylabel('Fairlead force [MN]','FontSize',15)
xlim([0 max(Time)]);
ax = gca;
ax.FontSize = 15;
legend({'l = 50 m','l = 25 m'},'FontSize',15,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\EndForceChain4TimeHistorySingleChainLowHigh_TEST1v2',NumLines);
print(str,'-dpng')

figure(34);
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
p1 = plot(Time,VonMisesTimeAllGeomtery(4,:)/1000,'-');
hold on 
p2 = plot(Time,VonMisesTimeAllGeomtery(1,:)/1000,'-');
grid minor
xlabel('Time [s]','FontSize',15)
ylabel('Fairlead stress [MPa]','FontSize',15)
xlim([0 max(Time)]);
ax = gca;
ax.FontSize = 15;
legend({'l = 50 m','l = 25 m'},'FontSize',15,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\EndForceChain4TimeHistorySingleVonMisesChainLowHigh_TEST1',NumLines);
print(str,'-dpng')

figure(57);
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
p1 = plot(Time,VonMisesTimeAllGeomtery(4,:)/1000,'--k');
hold on 
p2 = plot(Time,VonMisesTimeAllGeomtery(1,:)/1000,'-k');
grid minor
xlabel('Time [s]','FontSize',15)
ylabel('Fairlead stress [MPa]','FontSize',15)
xlim([0 max(Time)]);
ax = gca;
ax.FontSize = 15;
legend({'l = 50 m','l = 25 m'},'FontSize',15,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\EndForceChain4TimeHistorySingleVonMisesChainLowHigh_TEST1v2',NumLines);
print(str,'-dpng')

figure(35);
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
p1 = plot(Time,uTowerAllGeometry(4,:),'-');
hold on 
p2 = plot(Time,uTowerAllGeometry(1,:),'-');
grid minor
xlabel('Time [s]','FontSize',15)
ylabel('Hub displacement [m]','FontSize',15)
xlim([0 max(Time)]);
ax = gca;
ax.FontSize = 15;
legend({'l = 50 m','l = 25 m'},'FontSize',15,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\TowerDispTimeHistorySingleLowHigh_TEST1',NumLines);
print(str,'-dpng')


figure(55);
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
p1 = plot(Time,uTowerAllGeometry(4,:)*1000,'--k');
hold on 
p2 = plot(Time,uTowerAllGeometry(1,:)*1000,'-k');
grid minor
xlabel('Time [s]','FontSize',15)
ylabel('Hub displacement [mm]','FontSize',15)
xlim([0 max(Time)]);
ax = gca;
ax.FontSize = 15;
legend({'l = 50 m','l = 25 m'},'FontSize',15,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\TowerDispTimeHistorySingleLowHigh_TEST1v2',NumLines);
print(str,'-dpng')


figure(11);
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
p1 = plot(WallTensionMPAllGeometry(1).X,WallTensionMPAllGeometry(1).Max/1000,'-');
hold on 
p2 = plot(WallTensionMPAllGeometry(2).X,WallTensionMPAllGeometry(2).Max/1000,'-');
p3 = plot(WallTensionMPAllGeometry(3).X,WallTensionMPAllGeometry(3).Max/1000,'-');
p4 = plot(WallTensionMPAllGeometry(4).X,WallTensionMPAllGeometry(4).Max/1000,'-');
p5 = plot(WallTensionMPAllGeometry(5).X,WallTensionMPAllGeometry(5).Max/1000,'-');
p6 = plot(WallTensionMPAllGeometry(6).X,WallTensionMPAllGeometry(6).Max/1000,'-');
p7 = plot(WallTensionMPAllGeometry(7).X,WallTensionMPAllGeometry(7).Max/1000,'-');
p8 = plot(WallTensionMPAllGeometry(8).X,WallTensionMPAllGeometry(8).Max/1000,'-');
grid minor
xlim([0 TotalFoundationLength])
xlabel('L [m]','FontSize',15)
ylabel('Wall tension [MN]','FontSize',15)
ax = gca;
ax.FontSize = 15;
legend([p1 p2 p3 p4 p5 p6 p7 p8],{'l = 25 m','l = 30 m','l = 40 m','l = 50 m','l = 60 m','l = 70 m','l = 80 m','l = 90 m'},'FontSize',15,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_WallTension_TEST1',NumLines);
print(str,'-dpng')

figure(12);
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
p1 = plot(ExternalPressureMPAllGeometry(1).X,ExternalPressureMPAllGeometry(1).Max/1000,'-k');
hold on 
p2 = plot(ExternalPressureMPAllGeometry(2).X,ExternalPressureMPAllGeometry(2).Max/1000,'-k');
p3 = plot(ExternalPressureMPAllGeometry(3).X,ExternalPressureMPAllGeometry(3).Max/1000,'-k');
p4 = plot(ExternalPressureMPAllGeometry(4).X,ExternalPressureMPAllGeometry(4).Max/1000,'-k');
p5 = plot(ExternalPressureMPAllGeometry(5).X,ExternalPressureMPAllGeometry(5).Max/1000,'-k');
p6 = plot(ExternalPressureMPAllGeometry(6).X,ExternalPressureMPAllGeometry(6).Max/1000,'-k');
p7 = plot(ExternalPressureMPAllGeometry(7).X,ExternalPressureMPAllGeometry(7).Max/1000,'-k');
p8 = plot(ExternalPressureMPAllGeometry(8).X,ExternalPressureMPAllGeometry(8).Max/1000,'-k');
grid minor
xlim([0 TotalFoundationLength])
xlabel('L [m]','FontSize',15)
ylabel('External pressure [MPa]','FontSize',15)
ax = gca;
ax.FontSize = 15;
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_ExternalPressure_TEST1',NumLines);
print(str,'-dpng')


figure(13)
set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
y1 = [EndForceMaxAllGeometryChain4EndA(1,:)/1000];
x1 = [25];
y2 = [EndForceMaxAllGeometryChain4EndA(2,:)/1000];
x2 = [30];
y3 = [EndForceMaxAllGeometryChain4EndA(3,:)/1000];
x3 = [40];
y4 = [EndForceMaxAllGeometryChain4EndA(4,:)/1000];
x4 = [50];
y5 = [EndForceMaxAllGeometryChain4EndA(5,:)/1000];
x5 = [60];
y6 = [EndForceMaxAllGeometryChain4EndA(6,:)/1000];
x6 = [70];
y7 = [EndForceMaxAllGeometryChain4EndA(7,:)/1000];
x7 = [80];
y8 = [EndForceMaxAllGeometryChain4EndA(8,:)/1000];
x8 = [90];
bar(x1,y1,3.5)
hold on
bar(x2,y2,3.5)
bar(x3,y3,3.5)
bar(x4,y4,3.5)
bar(x5,y5,3.5)
bar(x6,y6,3.5)
bar(x7,y7,3.5)
bar(x8,y8,3.5)
bar(113,0)
grid minor
xlabel('l [m]','FontSize',15)
ylabel('Max chain force [MN]','FontSize',15)
ax = gca;
ax.FontSize = 15;
xticks([25 30 40 50 60 70 80 90 120])
xticklabels({'25','30','40','50','60','70','80','90',''})
legend({'L/l = 1.985','L/l = 1.745','L/l = 1.465','L/l = 1.316','L/l = 1.227','L/l = 1.170','L/l = 1.132','L/l = 1.105'},'FontSize',15,'Location','Northeast')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\Chain4_MaxEndForce_EndA_HingedBot_TEST1',NumLines);
print(str,'-dpng')



x0=14;
y0=2;

minXval = min([MomentMPAllGeometry(1).Max/1000,MomentMPAllGeometry(2).Max/1000,MomentMPAllGeometry(3).Max/1000,...
    MomentMPAllGeometry(4).Max/1000,MomentMPAllGeometry(5).Max/1000,MomentMPAllGeometry(6).Max/1000,...
    MomentMPAllGeometry(7).Max/1000,MomentMPAllGeometry(8).Max/1000,]);

maxXval = max([MomentMPAllGeometry(1).Max/1000,MomentMPAllGeometry(2).Max/1000,MomentMPAllGeometry(3).Max/1000,...
    MomentMPAllGeometry(4).Max/1000,MomentMPAllGeometry(5).Max/1000,MomentMPAllGeometry(6).Max/1000,...
    MomentMPAllGeometry(7).Max/1000,MomentMPAllGeometry(8).Max/1000,]);

figure(15);
set(gcf,'units','centimeters','position',[x0,y0,width,height])
p1 = plot(MomentMPAllGeometry(1).Max/1000,MomentMPAllGeometry(1).X,'-');
hold on 
p2 = plot(MomentMPAllGeometry(2).Max/1000,MomentMPAllGeometry(2).X,'-');
p3 = plot(MomentMPAllGeometry(3).Max/1000,MomentMPAllGeometry(3).X,'-');
p4 = plot(MomentMPAllGeometry(4).Max/1000,MomentMPAllGeometry(4).X,'-');
p5 = plot(MomentMPAllGeometry(5).Max/1000,MomentMPAllGeometry(5).X,'-');
p6 = plot(MomentMPAllGeometry(6).Max/1000,MomentMPAllGeometry(6).X,'-');
p7 = plot(MomentMPAllGeometry(7).Max/1000,MomentMPAllGeometry(7).X,'-');
p8 = plot(MomentMPAllGeometry(8).Max/1000,MomentMPAllGeometry(8).X,'-');

plot([minXval maxXval],[ch ch],':k'); % Chain connection dot line
str1 = {'Fairlead'};
text(6,ch+2,str1,'FontSize',16)

plot([minXval maxXval],[msl msl],':k'); % Chain connection dot line
str2 = {'MSL'};
text(6,msl+2,str2,'FontSize',16)

plot([minXval maxXval],[0 0],':k'); % Chain connection dot line
str3 = {'Seabed'};
text(6,0-2,str3,'FontSize',16)

plot([minXval maxXval],[yint yint],':k'); % Chain connection dot line
str3 = {'Interface'};
text(6,max(MomentMPAllGeometry(8).X)+2,str3,'FontSize',16)

grid minor
ylim([-5 75])
xlim([minXval maxXval])
ylabel('z [m]','FontSize',16)
xlabel('Moment [MNm]','FontSize',16)
ax = gca;
ax.FontSize = 16;
legend([p1 p2 p3 p4 p5 p6 p7 p8],{'l = 25 m','l = 30 m','l = 40 m','l = 50 m','l = 60 m','l = 70 m','l = 80 m','l = 90 m'},'FontSize',16,'Location','Southeast')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_BendingMoment_TEST1',NumLines);
print(str,'-dpng')

x0=28;
y0=2;

minXval = min([VonMisesMPAllGeometry(1).Max/1000,VonMisesMPAllGeometry(2).Max/1000,VonMisesMPAllGeometry(3).Max/1000,...
    VonMisesMPAllGeometry(4).Max/1000,VonMisesMPAllGeometry(5).Max/1000,VonMisesMPAllGeometry(6).Max/1000,...
    VonMisesMPAllGeometry(7).Max/1000,VonMisesMPAllGeometry(8).Max/1000,]);

maxXval = max([VonMisesMPAllGeometry(1).Max/1000,VonMisesMPAllGeometry(2).Max/1000,VonMisesMPAllGeometry(3).Max/1000,...
    VonMisesMPAllGeometry(4).Max/1000,VonMisesMPAllGeometry(5).Max/1000,VonMisesMPAllGeometry(6).Max/1000,...
    VonMisesMPAllGeometry(7).Max/1000,VonMisesMPAllGeometry(8).Max/1000,]);

figure(16);
set(gcf,'units','centimeters','position',[x0,y0,width,height])
p1 = plot(VonMisesMPAllGeometry(1).Max/1000,VonMisesMPAllGeometry(1).X,'-');
hold on 
p2 = plot(VonMisesMPAllGeometry(2).Max/1000,VonMisesMPAllGeometry(2).X,'-');
p3 = plot(VonMisesMPAllGeometry(3).Max/1000,VonMisesMPAllGeometry(3).X,'-');
p4 = plot(VonMisesMPAllGeometry(4).Max/1000,VonMisesMPAllGeometry(4).X,'-');
p5 = plot(VonMisesMPAllGeometry(5).Max/1000,VonMisesMPAllGeometry(5).X,'-');
p6 = plot(VonMisesMPAllGeometry(6).Max/1000,VonMisesMPAllGeometry(6).X,'-');
p7 = plot(VonMisesMPAllGeometry(7).Max/1000,VonMisesMPAllGeometry(7).X,'-');
p8 = plot(VonMisesMPAllGeometry(8).Max/1000,VonMisesMPAllGeometry(8).X,'-');

plot([minXval maxXval],[ch ch],':k'); % Chain connection dot line
str1 = {'Fairlead'};
text(29,ch+2,str1,'FontSize',16)

plot([minXval maxXval],[msl msl],':k'); % Chain connection dot line
str2 = {'MSL'};
text(29,msl+2,str2,'FontSize',16)

plot([minXval maxXval],[0 0],':k'); % Chain connection dot line
str3 = {'Seabed'};
text(32,0+2,str3,'FontSize',16)

plot([minXval maxXval],[yint yint],':k'); % Chain connection dot line
str3 = {'Interface'};
text(29,max(MomentMPAllGeometry(8).X)+2,str3,'FontSize',16)

grid minor
ylim([-5 75])
xlim([minXval maxXval])
ylabel('z [m]','FontSize',16)
xlabel('Max von Mises stress [MPa]','FontSize',16)
ax = gca;
ax.FontSize = 16;
legend([p1 p2 p3 p4 p5 p6 p7 p8],{'l = 25 m','l = 30 m','l = 40 m','l = 50 m','l = 60 m','l = 70 m','l = 80 m','l = 90 m'},'FontSize',16,'Location','Southeast')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_VonMises_TEST1',NumLines);
print(str,'-dpng')


minXval = min([ShearMPAllGeometry(1).Max/1000,ShearMPAllGeometry(2).Max/1000,ShearMPAllGeometry(3).Max/1000,...
    ShearMPAllGeometry(4).Max/1000,ShearMPAllGeometry(5).Max/1000,ShearMPAllGeometry(6).Max/1000,...
    ShearMPAllGeometry(7).Max/1000,ShearMPAllGeometry(8).Max/1000,]);

maxXval = max([ShearMPAllGeometry(1).Max/1000,ShearMPAllGeometry(2).Max/1000,ShearMPAllGeometry(3).Max/1000,...
    ShearMPAllGeometry(4).Max/1000,ShearMPAllGeometry(5).Max/1000,ShearMPAllGeometry(6).Max/1000,...
    ShearMPAllGeometry(7).Max/1000,ShearMPAllGeometry(8).Max/1000,]);

figure(17);
set(gcf,'units','centimeters','position',[x0,y0,width,height])
p1 = plot(ShearMPAllGeometry(1).Max/1000,ShearMPAllGeometry(1).X,'-');
hold on 
p2 = plot(ShearMPAllGeometry(2).Max/1000,ShearMPAllGeometry(2).X,'-');
p3 = plot(ShearMPAllGeometry(3).Max/1000,ShearMPAllGeometry(3).X,'-');
p4 = plot(ShearMPAllGeometry(4).Max/1000,ShearMPAllGeometry(4).X,'-');
p5 = plot(ShearMPAllGeometry(5).Max/1000,ShearMPAllGeometry(5).X,'-');
p6 = plot(ShearMPAllGeometry(6).Max/1000,ShearMPAllGeometry(6).X,'-');
p7 = plot(ShearMPAllGeometry(7).Max/1000,ShearMPAllGeometry(7).X,'-');
p8 = plot(ShearMPAllGeometry(8).Max/1000,ShearMPAllGeometry(8).X,'-');

plot([minXval maxXval],[ch ch],':k'); % Chain connection dot line
str1 = {'Fairlead'};
text(1.7,ch-2,str1,'FontSize',16)

plot([minXval maxXval],[msl msl],':k'); % Chain connection dot line
str2 = {'MSL'};
text(1.7,msl+2,str2,'FontSize',16)

plot([minXval maxXval],[0 0],':k'); % Chain connection dot line
str3 = {'Seabed'};
text(1.7,0+2,str3,'FontSize',16)

plot([minXval maxXval],[yint yint],':k'); % Chain connection dot line
str3 = {'Interface'};
text(1.7,max(MomentMPAllGeometry(8).X)+2,str3,'FontSize',16)

grid minor
ylim([-5 75])
xlim([minXval maxXval])
ylabel('z [m]','FontSize',16)
xlabel('Shear force [MN]','FontSize',16)
ax = gca;
ax.FontSize = 16;
legend([p1 p2 p3 p4 p5 p6 p7 p8],{'l = 25 m','l = 30 m','l = 40 m','l = 50 m','l = 60 m','l = 70 m','l = 80 m','l = 90 m'},'FontSize',16,'Location','Northeast')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_Shear_TEST1',NumLines);
print(str,'-dpng')

minXval = min([EffectiveTensionMPAllGeometry(1).Max/1000,EffectiveTensionMPAllGeometry(2).Max/1000,EffectiveTensionMPAllGeometry(3).Max/1000,...
    EffectiveTensionMPAllGeometry(4).Max/1000,EffectiveTensionMPAllGeometry(5).Max/1000,EffectiveTensionMPAllGeometry(6).Max/1000,...
    EffectiveTensionMPAllGeometry(7).Max/1000,EffectiveTensionMPAllGeometry(8).Max/1000,]);

maxXval = max([EffectiveTensionMPAllGeometry(1).Max/1000,EffectiveTensionMPAllGeometry(2).Max/1000,EffectiveTensionMPAllGeometry(3).Max/1000,...
    EffectiveTensionMPAllGeometry(4).Max/1000,EffectiveTensionMPAllGeometry(5).Max/1000,EffectiveTensionMPAllGeometry(6).Max/1000,...
    EffectiveTensionMPAllGeometry(7).Max/1000,EffectiveTensionMPAllGeometry(8).Max/1000,]);

figure(18);
set(gcf,'units','centimeters','position',[x0,y0,width,height])
p1 = plot(EffectiveTensionMPAllGeometry(1).Max/1000,EffectiveTensionMPAllGeometry(1).X,'-');
hold on 
p2 = plot(EffectiveTensionMPAllGeometry(2).Max/1000,EffectiveTensionMPAllGeometry(2).X,'-');
p3 = plot(EffectiveTensionMPAllGeometry(3).Max/1000,EffectiveTensionMPAllGeometry(3).X,'-');
p4 = plot(EffectiveTensionMPAllGeometry(4).Max/1000,EffectiveTensionMPAllGeometry(4).X,'-');
p5 = plot(EffectiveTensionMPAllGeometry(5).Max/1000,EffectiveTensionMPAllGeometry(5).X,'-');
p6 = plot(EffectiveTensionMPAllGeometry(6).Max/1000,EffectiveTensionMPAllGeometry(6).X,'-');
p7 = plot(EffectiveTensionMPAllGeometry(7).Max/1000,EffectiveTensionMPAllGeometry(7).X,'-');
p8 = plot(EffectiveTensionMPAllGeometry(8).Max/1000,EffectiveTensionMPAllGeometry(8).X,'-');

plot([minXval maxXval],[ch ch],':k'); % Chain connection dot line
str1 = {'Fairlead'};
text(-28,ch-2,str1,'FontSize',16)

plot([minXval maxXval],[msl msl],':k'); % Chain connection dot line
str2 = {'MSL'};
text(-28,msl+2,str2,'FontSize',16)

plot([minXval maxXval],[0 0],':k'); % Chain connection dot line
str3 = {'Seabed'};
text(-28,0+2,str3,'FontSize',16)

plot([minXval maxXval],[yint yint],':k'); % Chain connection dot line
str3 = {'Interface'};
text(-28,max(EffectiveTensionMPAllGeometry(8).X)+2,str3,'FontSize',16)

grid minor
%title(sprintf('MP reaction force at end A. Hinged bottom, static load angle %d deg',LoadAngle))
ylim([-5 75])
xlim([minXval maxXval])
ylabel('z [m]','FontSize',16)
xlabel('Effective tension [MN]','FontSize',16)
ax = gca;
ax.FontSize = 16;
legend([p1 p2 p3 p4 p5 p6 p7 p8],{'l = 25 m','l = 30 m','l = 40 m','l = 50 m','l = 60 m','l = 70 m','l = 80 m','l = 90 m'},'FontSize',16,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_EffectiveTension_TEST1',NumLines);
print(str,'-dpng')


minXval = min([WallTensionMPAllGeometry(1).Max/1000,WallTensionMPAllGeometry(2).Max/1000,WallTensionMPAllGeometry(3).Max/1000,...
    WallTensionMPAllGeometry(4).Max/1000,WallTensionMPAllGeometry(5).Max/1000,WallTensionMPAllGeometry(6).Max/1000,...
    WallTensionMPAllGeometry(7).Max/1000,WallTensionMPAllGeometry(8).Max/1000,]);

maxXval = max([WallTensionMPAllGeometry(1).Max/1000,WallTensionMPAllGeometry(2).Max/1000,WallTensionMPAllGeometry(3).Max/1000,...
    WallTensionMPAllGeometry(4).Max/1000,WallTensionMPAllGeometry(5).Max/1000,WallTensionMPAllGeometry(6).Max/1000,...
    WallTensionMPAllGeometry(7).Max/1000,WallTensionMPAllGeometry(8).Max/1000,]);

figure(19);
set(gcf,'units','centimeters','position',[x0,y0,width,height])
p1 = plot(WallTensionMPAllGeometry(1).Max/1000,WallTensionMPAllGeometry(1).X,'-');
hold on 
p2 = plot(WallTensionMPAllGeometry(2).Max/1000,WallTensionMPAllGeometry(2).X,'-');
p3 = plot(WallTensionMPAllGeometry(3).Max/1000,WallTensionMPAllGeometry(3).X,'-');
p4 = plot(WallTensionMPAllGeometry(4).Max/1000,WallTensionMPAllGeometry(4).X,'-');
p5 = plot(WallTensionMPAllGeometry(5).Max/1000,WallTensionMPAllGeometry(5).X,'-');
p6 = plot(WallTensionMPAllGeometry(6).Max/1000,WallTensionMPAllGeometry(6).X,'-');
p7 = plot(WallTensionMPAllGeometry(7).Max/1000,WallTensionMPAllGeometry(7).X,'-');
p8 = plot(WallTensionMPAllGeometry(8).Max/1000,WallTensionMPAllGeometry(8).X,'-');

plot([minXval maxXval],[ch ch],':k'); % Chain connection dot line
str1 = {'Fairlead'};
text(-33,ch-2,str1,'FontSize',16)

plot([minXval maxXval],[msl msl],':k'); % Chain connection dot line
str2 = {'MSL'};
text(-33,msl+2,str2,'FontSize',16)

plot([minXval maxXval],[0 0],':k'); % Chain connection dot line
str3 = {'Seabed'};
text(-33,0+2,str3,'FontSize',16)

plot([minXval maxXval],[yint yint],':k'); % Chain connection dot line
str3 = {'Interface'};
text(-33,max(WallTensionMPAllGeometry(8).X)+2,str3,'FontSize',16)

grid minor
ylim([-5 75])
xlim([minXval maxXval])
ylabel('z [m]','FontSize',16)
xlabel('Wall tension [MN]','FontSize',16)
ax = gca;
ax.FontSize = 16;
legend([p1 p2 p3 p4 p5 p6 p7 p8],{'l = 25 m','l = 30 m','l = 40 m','l = 50 m','l = 60 m','l = 70 m','l = 80 m','l = 90 m'},'FontSize',16,'Location','Northwest')
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_WallTension_TEST1',NumLines);
print(str,'-dpng')


minXval = min([ExternalPressureMPAllGeometry(1).Max/1000,ExternalPressureMPAllGeometry(2).Max/1000,ExternalPressureMPAllGeometry(3).Max/1000,...
    ExternalPressureMPAllGeometry(4).Max/1000,ExternalPressureMPAllGeometry(5).Max/1000,ExternalPressureMPAllGeometry(6).Max/1000,...
    ExternalPressureMPAllGeometry(7).Max/1000,ExternalPressureMPAllGeometry(8).Max/1000,]);

maxXval = max([ExternalPressureMPAllGeometry(1).Max/1000,ExternalPressureMPAllGeometry(2).Max/1000,ExternalPressureMPAllGeometry(3).Max/1000,...
    ExternalPressureMPAllGeometry(4).Max/1000,ExternalPressureMPAllGeometry(5).Max/1000,ExternalPressureMPAllGeometry(6).Max/1000,...
    ExternalPressureMPAllGeometry(7).Max/1000,ExternalPressureMPAllGeometry(8).Max/1000,]);

figure(21);
set(gcf,'units','centimeters','position',[x0,y0,width,height])
p1 = plot(ExternalPressureMPAllGeometry(1).Max/1000,ExternalPressureMPAllGeometry(1).X,'-k');
hold on 
p2 = plot(ExternalPressureMPAllGeometry(2).Max/1000,ExternalPressureMPAllGeometry(2).X,'-k');
p3 = plot(ExternalPressureMPAllGeometry(3).Max/1000,ExternalPressureMPAllGeometry(3).X,'-k');
p4 = plot(ExternalPressureMPAllGeometry(4).Max/1000,ExternalPressureMPAllGeometry(4).X,'-k');
p5 = plot(ExternalPressureMPAllGeometry(5).Max/1000,ExternalPressureMPAllGeometry(5).X,'-k');
p6 = plot(ExternalPressureMPAllGeometry(6).Max/1000,ExternalPressureMPAllGeometry(6).X,'-k');
p7 = plot(ExternalPressureMPAllGeometry(7).Max/1000,ExternalPressureMPAllGeometry(7).X,'-k');
p8 = plot(ExternalPressureMPAllGeometry(8).Max/1000,ExternalPressureMPAllGeometry(8).X,'-k');

plot([minXval maxXval],[ch ch],':k'); % Chain connection dot line
str1 = {'Fairlead'};
text(0.35,ch+2,str1,'FontSize',16)

plot([minXval maxXval],[msl msl],':k'); % Chain connection dot line
str2 = {'MSL'};
text(0.35,msl+2,str2,'FontSize',16)

plot([minXval maxXval],[0 0],':k'); % Chain connection dot line
str3 = {'Seabed'};
text(0.35,0+2,str3,'FontSize',16)

plot([minXval maxXval],[yint yint],':k'); % Chain connection dot line
str3 = {'Interface'};
text(0.35,max(ExternalPressureMPAllGeometry(8).X)+2,str3,'FontSize',16)

grid minor
ylim([-5 75])
xlim([minXval maxXval])
ylabel('z [m]','FontSize',16)
xlabel('External pressure [MPa]','FontSize',16)
ax = gca;
ax.FontSize = 16;
str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_ExternalPressure_TEST1',NumLines);
print(str,'-dpng')


%%%% THESE FIGURES MAY BE DELETED %%%%%%%

% figure(31);
% set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
% p1 = plot(Time,EndForceTimeHistory(1,:)/1000,'-');
% hold on 
% grid minor
% xlabel('Time [s]','FontSize',16)
% ylabel('Fairlead force [MN]','FontSize',16)
% xlim([0 max(Time)]);
% ax = gca;
% ax.FontSize = 16;
% legend({'L/l = 1.985'},'FontSize',16,'Location','Northwest')
% str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\EndForceChain4TimeHistorySingleChainHigh_TEST1',NumLines);
% print(str,'-dpng')

% figure(32);
% set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
% p1 = plot(Time,EndForceTimeHistory(4,:)/1000,'-');
% hold on 
% grid minor
% xlabel('Time [s]','FontSize',16)
% ylabel('Fairlead force [MN]','FontSize',16)
% xlim([0 max(Time)]);
% ax = gca;
% ax.FontSize = 16;
% legend({'L/l = 1.316'},'FontSize',16,'Location','Northwest')
% str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\EndForceChain4TimeHistorySingleChainLow_TEST1',NumLines);
% print(str,'-dpng')


% figure(4);
% set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
% p1 = plot(MomentMPAllGeometry(1).X,MomentMPAllGeometry(1).Max/1000,'-');
% hold on 
% p2 = plot(MomentMPAllGeometry(2).X,MomentMPAllGeometry(2).Max/1000,'-');
% p3 = plot(MomentMPAllGeometry(3).X,MomentMPAllGeometry(3).Max/1000,'-');
% p4 = plot(MomentMPAllGeometry(4).X,MomentMPAllGeometry(4).Max/1000,'-');
% p5 = plot(MomentMPAllGeometry(5).X,MomentMPAllGeometry(5).Max/1000,'-');
% p6 = plot(MomentMPAllGeometry(6).X,MomentMPAllGeometry(6).Max/1000,'-');
% p7 = plot(MomentMPAllGeometry(7).X,MomentMPAllGeometry(7).Max/1000,'-');
% p8 = plot(MomentMPAllGeometry(8).X,MomentMPAllGeometry(8).Max/1000,'-');
% grid minor
% xlim([0 TotalFoundationLength])
% xlabel('L [m]','FontSize',16)
% ylabel('Moment [MNm]','FontSize',16)
% ax = gca;
% ax.FontSize = 16;
% legend([p1 p2 p3 p4 p5 p6 p7 p8],{'l = 25 m','l = 30 m','l = 40 m','l = 50 m','l = 60 m','l = 70 m','l = 80 m','l = 90 m'},'FontSize',16,'Location','Northwest')
% str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_BendingMoment_TEST1',NumLines);
% print(str,'-dpng')



% figure(5);
% set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
% p1 = plot(VonMisesMPAllGeometry(1).X,VonMisesMPAllGeometry(1).Max/1000,'-');
% hold on 
% p2 = plot(VonMisesMPAllGeometry(2).X,VonMisesMPAllGeometry(2).Max/1000,'-');
% p3 = plot(VonMisesMPAllGeometry(3).X,VonMisesMPAllGeometry(3).Max/1000,'-');
% p4 = plot(VonMisesMPAllGeometry(4).X,VonMisesMPAllGeometry(4).Max/1000,'-');
% p5 = plot(VonMisesMPAllGeometry(5).X,VonMisesMPAllGeometry(5).Max/1000,'-');
% p6 = plot(VonMisesMPAllGeometry(6).X,VonMisesMPAllGeometry(6).Max/1000,'-');
% p7 = plot(VonMisesMPAllGeometry(7).X,VonMisesMPAllGeometry(7).Max/1000,'-');
% p8 = plot(VonMisesMPAllGeometry(8).X,VonMisesMPAllGeometry(8).Max/1000,'-');
% grid minor
% xlim([0 TotalFoundationLength])
% xlabel('L [m]','FontSize',16)
% ylabel('Max von Mises stress [MPa]','FontSize',16)
% ax = gca;
% ax.FontSize = 16;
% legend([p1 p2 p3 p4 p5 p6 p7 p8],{'l = 25 m','l = 30 m','l = 40 m','l = 50 m','l = 60 m','l = 70 m','l = 80 m','l = 90 m'},'FontSize',16,'Location','Northwest')
% str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_VonMises_TEST1',NumLines);
% print(str,'-dpng')




% figure(6);
% set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
% p1 = plot(ShearMPAllGeometry(1).X,ShearMPAllGeometry(1).Max/1000,'-');
% hold on 
% p2 = plot(ShearMPAllGeometry(2).X,ShearMPAllGeometry(2).Max/1000,'-');
% p3 = plot(ShearMPAllGeometry(3).X,ShearMPAllGeometry(3).Max/1000,'-');
% p4 = plot(ShearMPAllGeometry(4).X,ShearMPAllGeometry(4).Max/1000,'-');
% p5 = plot(ShearMPAllGeometry(5).X,ShearMPAllGeometry(5).Max/1000,'-');
% p6 = plot(ShearMPAllGeometry(6).X,ShearMPAllGeometry(6).Max/1000,'-');
% p7 = plot(ShearMPAllGeometry(7).X,ShearMPAllGeometry(7).Max/1000,'-');
% p8 = plot(ShearMPAllGeometry(8).X,ShearMPAllGeometry(8).Max/1000,'-');
% grid minor
% xlim([0 TotalFoundationLength])
% xlabel('L [m]','FontSize',16)
% ylabel('Shear force [MN]','FontSize',16)
% ax = gca;
% ax.FontSize = 16;
% legend([p1 p2 p3 p4 p5 p6 p7 p8],{'l = 25 m','l = 30 m','l = 40 m','l = 50 m','l = 60 m','l = 70 m','l = 80 m','l = 90 m'},'FontSize',16,'Location','Northeast')
% str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_Shear_TEST1',NumLines);
% print(str,'-dpng')




% figure(9)
% set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
% y = [6.348 1.328; 6.127 1.381; 5.813 1.290; 5.730 1.399; 5.695 1.424; 5.702 1.474; 5.759 1.659; 5.821 1.798];
% x = [25 30 40 50 60 70 80 90];
% b = bar(x,y,0.9);
% hold on
% b(2).FaceColor = 'yellow';
% bar(103,0)
% grid minor
% xlabel('l [m]','FontSize',16)
% ylabel('Period [s]','FontSize',16)
% ax = gca;
% ax.FontSize = 16;
% xticks([25 30 40 50 60 70 80 90])
% xticklabels({'25','30','40','50','60','70','80','90'})
% legend({'First period','Second period'},'FontSize',16,'Location','Northeast')
% str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\Periods_TEST1',NumLines);
% print(str,'-dpng')



% figure(10);
% set(gcf,'units','centimeters','position',[x0,y0,widthNEW,heightNEW])
% p1 = plot(EffectiveTensionMPAllGeometry(1).X,EffectiveTensionMPAllGeometry(1).Max/1000,'-');
% hold on 
% p2 = plot(EffectiveTensionMPAllGeometry(2).X,EffectiveTensionMPAllGeometry(2).Max/1000,'-');
% p3 = plot(EffectiveTensionMPAllGeometry(3).X,EffectiveTensionMPAllGeometry(3).Max/1000,'-');
% p4 = plot(EffectiveTensionMPAllGeometry(4).X,EffectiveTensionMPAllGeometry(4).Max/1000,'-');
% p5 = plot(EffectiveTensionMPAllGeometry(5).X,EffectiveTensionMPAllGeometry(5).Max/1000,'-');
% p6 = plot(EffectiveTensionMPAllGeometry(6).X,EffectiveTensionMPAllGeometry(6).Max/1000,'-');
% p7 = plot(EffectiveTensionMPAllGeometry(7).X,EffectiveTensionMPAllGeometry(7).Max/1000,'-');
% p8 = plot(EffectiveTensionMPAllGeometry(8).X,EffectiveTensionMPAllGeometry(8).Max/1000,'-');
% grid minor
% xlim([0 TotalFoundationLength])
% xlabel('L [m]','FontSize',16)
% ylabel('Effective tension [MN]','FontSize',16)
% ax = gca;
% ax.FontSize = 16;
% legend([p1 p2 p3 p4 p5 p6 p7 p8],{'l = 25 m','l = 30 m','l = 40 m','l = 50 m','l = 60 m','l = 70 m','l = 80 m','l = 90 m'},'FontSize',16,'Location','Northwest')
% str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\DynamicAnalysisOfBuoyancyGuyedMonopiles\\FIGURES\\%sLines\\JONSWAP\\MP_EffectiveTension_TEST1',NumLines);
% print(str,'-dpng')








%%% To find a specific values %%%
% Find the time two which the highest water level occurs

%k = find(abs(etaAllGeometry(1,:)-max(etaAllGeometry(1,:))) < 0.00001)
%TimeStepSize = 0.01;
%TimeOfMaxWave = k*TimeStepSize






