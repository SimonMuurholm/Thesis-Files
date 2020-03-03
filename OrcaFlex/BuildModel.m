function [NumLines,MSL,LoadAngle,lvals,Lvals,Time,umaxMP,l,L,EndForceResults,uMP,... 
    EndForce,umaxMPAllGeometry,EndForceMaxAllGeometryChain1,EndForceMaxAllGeometryMP,... 
    MomentMPAllGeometry,VonMisesMPAllGeometry,ShearMPAllGeometry,umaxTowerAllGeometry,... 
    EndForceMaxAllGeometryChain4,etaAllGeometry,EffectiveTensionMPAllGeometry,...
    WallTensionMPAllGeometry,ExternalPressureMPAllGeometry,TotalFoundationLength,...
    EndForceMaxAllGeometryChain4EndA,EndForceTimeHistory,VonMisesTimeAllGeomtery,...
    uTowerAllGeometry] = BuildModel(filename)

NumGeomInputFiles = 8; % Number of footprint lengths i.e input files
NumModels = 1; % Number of models run from number of sheets in geometry input 
EndStageLength = 600; % Set length of simulation
MSL = 50; % Mean sea level
BGM = 46;
TPLength = 15.7;
SkirtLength = 10.5;
TotalFoundationLength = BGM + TPLength + SkirtLength;
LoadAngle = 0; % Angle of static load, for three lines model 0 or 60; For six lines model 0 or 30
ChainSegmentLength = 2;
StaticLoad = 0; % 1 for applied load, 0 for no applied load
NumLines = 'TOWERSix'; % Number of lines in BGM design ie. 'Six' or 'Three'
TowerInput = '13MW_GE_TowerGeometry.xlsx'; % File containing tower geometry
TowerGeometry = readtable(TowerInput,'sheet','Sheet1','ReadRowNames',false); % Tower
MPInput = 'MP_Geometry.xlsx'; % File containing MP geometry
MPGeometry = readtable(MPInput,'sheet','Sheet1','ReadRowNames',false); % MP
TowerMassesInput = 'TowerMassesv2.xlsx'; % File containing MP geometry
TowerMasses = readtable(TowerMassesInput,'sheet','Sheet1','ReadRowNames',false); % MP
DampingInput = 'Damping.xlsx'; % File containing MP geometry
Damping = readtable(DampingInput,'sheet','Sheet1','ReadRowNames',false); % MP
for jj = 1:NumGeomInputFiles

    strLines = sprintf('%sLineModelToLoad.dat',NumLines);
    model = ofxModel; % Create an empty model
    model.LoadData(strLines) % Load in geometry from pre-made file
    

                        %%%% INPUT FILES HERE %%%%
    Wavetype = 'JONSWAP';
       switch(Wavetype)
   case 'Airy'
       WaveInput = 'WaveInputFileAiry.xlsx'; % Regular waves
       case 'JONSWAP'
        WaveInput = 'WaveInputFileJONSWAP.xlsx'; % Irregular waves "WaveInputFileJONSWAP_STATIC" for no waves
       end
      
    str = sprintf('D:\\Dropbox (WoodThilsted)\\!E - Employees drive\\SMH\\Master Thesis\\ScriptedOrcaFlex\\%sLineModel\\InputFilesJONSWAP\\GeometryInputFileHINGE_JONSWAP',NumLines);
    baseGeomInputFileName = str; % Takes in input files from path above
    GeomInput = strcat(baseGeomInputFileName, sprintf('%d',jj),'.xlsx'); 
    WaveEnv = readtable(WaveInput,'sheet','Waves','ReadRowNames',false); % Envelope containing loads
    Lines = readtable(GeomInput,'sheet','Lines_Sheet','ReadRowNames',false); % Envelope containing loads 
    LineTypeName = readtable(GeomInput,'sheet','LineTypeName_Sheet','ReadRowNames',false); 
    Buoys = readtable(GeomInput,'sheet','Buoys_Sheet','ReadRowNames',false); % Buoys

    
    % Dynamics
    LoggingStepSize = 0.01;
    TimeStepSize = 0.01;
    model.general.TargetLogSampleInterval = LoggingStepSize;
    model.general.ImplicitConstantTimeStep = TimeStepSize;
    RayleighDamping = model('Damping1');
    
    % Stages 

    model.general.StageDuration(1) = 0.01; % indexed data, note stage 1 is the Build-up
    model.general.StageDuration(2) = EndStageLength;

                        %%%% VARIABLE DATA %%%%
                        
    % Box load  
    
    BoxX = model.CreateObject(ofx.otLoadForce, 'BoxX');
    BoxX.NumberOfRows = 3;
    BoxX.IndependentValue(1) = 0;
    BoxX.DependentValue(1) = 1000*cos(LoadAngle*(pi/180))*StaticLoad;
    BoxX.IndependentValue(2) = 0.5;
    BoxX.DependentValue(2) = 1000*cos(LoadAngle*(pi/180))*StaticLoad;
    BoxX.IndependentValue(3) = 1;
    BoxX.DependentValue(3) = 1000*cos(LoadAngle*(pi/180))*StaticLoad;
    
    BoxY = model.CreateObject(ofx.otLoadForce, 'BoxY');
    BoxY.NumberOfRows = 3;
    BoxY.IndependentValue(1) = 0;
    BoxY.DependentValue(1) = 1000*sin(LoadAngle*(pi/180))*StaticLoad;
    BoxY.IndependentValue(2) = 0.5;
    BoxY.DependentValue(2) = 1000*sin(LoadAngle*(pi/180))*StaticLoad;
    BoxY.IndependentValue(3) = 1;
    BoxY.DependentValue(3) = 1000*sin(LoadAngle*(pi/180))*StaticLoad;
    
    % Zero load                    
    Zero = model.CreateObject(ofx.otLoadForce, 'Zero');
    Zero.NumberOfRows = 1;
    Zero.IndependentValue(1) = 0;
    Zero.DependentValue(1) = 0;

    % Linetypes outside loop

for c = 2:size(LineTypeName,1) % Create all linetypes needed for chains. Note loop starts in 2 because linetype1 exists in input Orca file
    linetype = model.CreateObject(ofx.otLineType, sprintf('linetype%d',c))
end

for b = 1:size(TowerGeometry,1) % Create all linetypes needed for tower segmentation
    linetype = model.CreateObject(ofx.otLineType, sprintf('Tower%d',b))
end

for d = 1:size(MPGeometry,1) % Create all linetypes needed for MP segmentation
    linetype = model.CreateObject(ofx.otLineType, sprintf('MP%d',d))
end
 
    numArrays = 7; % One row for each line in model. Multiplum of 4 if more than disp wanted
    DisplacementResults = cell(numArrays,NumModels); % Each column is another model
    EndForceResults = cell(numArrays,NumModels); % Each column is another model
    MomentResults = cell(numArrays,NumModels);
    VonMisesResults = cell(numArrays,NumModels);
    ShearResults = cell(numArrays,NumModels);
    DisplacementTowerResults = cell(numArrays,NumModels);
    WaveElevationResults = cell(numArrays,NumModels);
    VonMisesResultsTime = cell(numArrays,NumModels);
    EffectiveTensionResults = cell(numArrays,NumModels);
    WallTensionResults = cell(numArrays,NumModels);
    ExternalPressureResults = cell(numArrays,NumModels);
    WholeSystemFrequencyResults1 = cell(1,NumModels); % Each column is another model, 1 because frequency found for whole system
    WholeSystemFrequencyResults2 = cell(1,NumModels);
    WholeSystemFrequencyResults3 = cell(1,NumModels);
for ii = 1:NumModels % Geometry loop over number of sheets
    baseSheetName = 'LinesBC_Sheet';
    extension = num2str(ii);
    SheetName = strcat(baseSheetName, extension);
    LinesBC = readtable(GeomInput,'sheet',SheetName,'ReadRowNames',false);
    
                            %%%% ENVIRONMENT %%%%
    
    model.environment.WaterDepth = MSL; % [m]
    
    % Sea density
    model.environment.Density = 1; % [t/m^3]
    
                        
    % Waves
    environment = model.environment;
    environment.NumberOfWaveTrains = 1;
   
    environment.WaveName(1) = cell2mat(WaveEnv.WaveName(1));
    environment.WaveTrainType = cell2mat(WaveEnv.WaveTrainType(1));
    environment.WaveTrainDirection = WaveEnv.WaveTrainDirection(1);
    
    environment.WaveJONSWAPParameters = 'Automatic';
     
   switch(Wavetype)
   case 'Airy'
           % Airy (Regular)
    environment.WaveHeight = WaveEnv.WaveHeight(1);
    environment.WavePeriod = WaveEnv.WavePeriod(1);
       
       case 'JONSWAP'
               % JONSWAP (Irregular)
    environment.WaveTrainHs = WaveEnv.WaveTrainHs(1);                  
    environment.WaveTrainTz = WaveEnv.WaveTrainTz(1);
   end

    
    for i = 1:size(Lines,1) % Loop over number of lines to build model
        
                        %%%% LINE DATA %%%% 
           
           line = model(cell2mat(Lines{i,'Lines'}))
           linetype = model(cell2mat(LineTypeName{i,'Line_Types'}))      
           linetype.Category = cell2mat(LineTypeName.Category(i));
           linetype.WizardCalculation = cell2mat(LineTypeName.WizardCalculation(i));
           linetype.ChainBarDiameter = LineTypeName.ChainBarDiameter(i);
           linetype.PipeOuterDiameter = LineTypeName.PipeOuterDiameter(i);
           linetype.PipeWallThickness = LineTypeName.PipeWallThickness(i);
           
           RayleighDamping.DampingRatio = Damping.DampingRatio(jj);
           RayleighDamping.Period1 = Damping.Period1(jj);
           RayleighDamping.Period2 = Damping.Period2(jj);
           linetype.RayleighDampingCoefficients = 'Damping1';
           
           linetype.InvokeWizard;
    
    % Set the Line end positions
           line.EndAX = LinesBC.EndAX(i);
           line.EndAY = LinesBC.EndAY(i);
           line.EndAZ = LinesBC.EndAZ(i);
           line.EndBZ = LinesBC.EndBZ(i);
              
    % Set the line end connections stiffness
           line.EndAxBendingStiffness = LinesBC.EndAxBendingStiffness(i);
           line.EndAyBendingStiffness = LinesBC.EndAyBendingStiffness(i);
    % Give the line sections
           line.NumberOfSections = 1; 
    % Set the section lengths and segmentation, first of all section 1 ...
           line.Length(1) = LinesBC.Length(i);
%            line.TargetSegmentLength(1) = LinesBC.TargetSegmentLength(i);
           line.TargetSegmentLength(1) = ChainSegmentLength;
           line.LineType(1) = cell2mat(LineTypeName.LineType(i));
     % Exclude seabed friction      
           line.IncludeSeabedFrictionInStatics = cell2mat(LinesBC.IncludeSeabedFrictionInStatics(i));
           
           
                % Nacelle (attached mass)
       if i == 1
           buoy = model(cell2mat(Buoys{i,'Buoys_3D'}))
           buoy.InitialZ = LinesBC.InitialZ(i);
           buoy.Mass = LinesBC.Mass(i);          
       else
           continue
       end
           
     % Tower masses
       if i == 1
           
           for g = 1:size(TowerMasses,1) % add more buoys (masses) to MP line (tower)
                buoy = model(cell2mat(TowerMasses{g,'Buoys_3D'}))         
                buoy.InitialZ = TowerMasses.InitialZ(g);
                buoy.Mass = TowerMasses.Mass(g); 
           end
       else
           continue
       end
       
       
       %%%% MP %%%%
       
           if i == 1                         
           for a = 1:size(MPGeometry,1) % add more linetypes to MP line
                line.NumberOfSections = size(MPGeometry,1);
                line = model(cell2mat(Lines{i,'Lines'})) % Reads in MP every time
                linetype = model(cell2mat(MPGeometry{a,'Line_Types'}));
                line.LineType(a) = cell2mat(MPGeometry.LineType(a));
                linetype.Category = cell2mat(MPGeometry.Category(a));
                linetype.WizardCalculation = cell2mat(MPGeometry.WizardCalculation(a));
                linetype.PipeOuterDiameter = MPGeometry.PipeOuterDiameter(a);
                linetype.PipeWallThickness = MPGeometry.PipeWallThickness(a);
                
                RayleighDamping.DampingRatio = Damping.DampingRatio(jj);
                RayleighDamping.Period1 = Damping.Period1(jj);
                RayleighDamping.Period2 = Damping.Period2(jj);
                linetype.RayleighDampingCoefficients = 'Damping1';
                
                linetype.InvokeWizard;  
                line.Length(a) = MPGeometry.SegmentLength(a);
                line.TargetSegmentLength(a) = MPGeometry.TargetSegmentLength(a);
                
           end
           
       else
           continue
       end     
       
           
       %%%%% TOWER %%%%
       
       if i == 1          
           for a = 1:size(TowerGeometry,1) % add more linetypes to MP line
                line.NumberOfSections = size(TowerGeometry,1)+size(MPGeometry,1);
                line = model(cell2mat(Lines{i,'Lines'})) % Reads in MP every time
                linetype = model(cell2mat(TowerGeometry{a,'Line_Types'}));
                line.LineType(a+size(MPGeometry,1)) = cell2mat(TowerGeometry.LineType(a)); % Assign linetype 
                linetype.Category = cell2mat(TowerGeometry.Category(a));
                linetype.WizardCalculation = cell2mat(TowerGeometry.WizardCalculation(a));
                linetype.PipeOuterDiameter = TowerGeometry.PipeOuterDiameter(a);
                linetype.PipeWallThickness = TowerGeometry.PipeWallThickness(a);
                
                RayleighDamping.DampingRatio = Damping.DampingRatio(jj);
                RayleighDamping.Period1 = Damping.Period1(jj);
                RayleighDamping.Period2 = Damping.Period2(jj);
                linetype.RayleighDampingCoefficients = 'Damping1';
                
                linetype.InvokeWizard;   
                line.Length(a+size(MPGeometry,1)) = TowerGeometry.SegmentLength(a);
                line.TargetSegmentLength(a+size(MPGeometry,1)) = TowerGeometry.TargetSegmentLength(a); % For descrization of line segment
           end
           
       else
           continue
       end
       

                         %     % Add a variable load
    LoadStart = 0;
    LoadEnd = EndStageLength;   
    Time = [LoadStart:LoggingStepSize:LoadEnd];
    
    line.NumberOfAppliedLoads = 1;
    line.AppliedLoadOriginZ(1) = MSL;
    line.AppliedLoadzRel(1) = 'End A';
    line.TopEnd = 'End B';
    line.AppliedForceX(1) = cell2mat(LinesBC.AppliedForceX(i));
    line.AppliedForceY(1) = cell2mat(LinesBC.AppliedForceY(i));

    end

                            %%%% CALCULATE STATICS AND MODESHAPES %%%%

%      model.CalculateStatics; % Run static simulation for modal analysis
%      modes = ofxModes(line, ofxModes.ModalAnalysisSpecification(true, -1, -1, true)); % Calculate all modes
%      mode1 = modes.modeDetails(1); % Get details of first mode
%      period1 = mode1.period % a single floating point value
%      mode2 = modes.modeDetails(3); % Get details of first mode
%      period2 = mode2.period % a single floating point value
%      mode3 = modes.modeDetails(5); % Get details of first mode
%      period3 = mode3.period % a single floating point value
     
%     shapeWrtGlobal = mode.shapeWrtGlobal; % an array of length dofCount
%     modeType = mode.modeType; % one of the pre-defined values ofx.mtTransverse, ofx.mtInline, etc.
%   
    
                    %%%% SAVING DATA %%%%
    fprintf('Running simulation \n')  
    model.RunSimulation; % Run dynamic simulation
    model.SaveSimulation(strcat(filename,sprintf('Geometry%d',ii),'.sim'));
    
    model.SaveData(strcat(filename,sprintf('Geometry%d',ii),'.dat'));
    fprintf('Geometry loop finished \n')
    model.LoadSimulation(strcat(filename,sprintf('Geometry%d',ii),'.sim'))

    
                    %%%% PRODUCE RESULTS %%%%

    MPline = model('MP');
    Chain1line = model('Chain1');
    Chain2line = model('Chain2');
    Chain3line = model('Chain3');
    Chain4line = model('Chain4');
%     Chain5line = model('Chain5');
%     Chain6line = model('Chain6');
    
    % The following lines are placed into cells, ii refers to geometry file
    DisplacementResults{1,ii} = MPline.TimeHistory('X',ofx.Period(0, EndStageLength), ofx.oeArcLength(MSL));
    DisplacementResults{2,ii} = Chain1line.TimeHistory('X',ofx.Period(0, EndStageLength), ofx.oeEndB);
    DisplacementResults{3,ii} = Chain2line.TimeHistory('X',ofx.Period(0, EndStageLength), ofx.oeEndB);
    DisplacementResults{4,ii} = Chain3line.TimeHistory('X',ofx.Period(0, EndStageLength), ofx.oeEndB);
%     DisplacementResults{5,ii} = Chain4line.TimeHistory('X',ofx.Period(0, EndStageLength), ofx.oeEndB);
%     DisplacementResults{6,ii} = Chain5line.TimeHistory('X',ofx.Period(0, EndStageLength), ofx.oeEndB);
%     DisplacementResults{7,ii} = Chain6line.TimeHistory('X',ofx.Period(0, EndStageLength), ofx.oeEndB);


    EndForceResults{1,ii} = MPline.TimeHistory('End GZ Force',ofx.Period(0, EndStageLength), ofx.oeEndA);
    EndForceResults{2,ii} = Chain1line.TimeHistory('End Force',ofx.Period(0, EndStageLength), ofx.oeEndB);
    EndForceResults{3,ii} = Chain2line.TimeHistory('End Force',ofx.Period(0, EndStageLength), ofx.oeEndB);
    EndForceResults{4,ii} = Chain3line.TimeHistory('End Force',ofx.Period(0, EndStageLength), ofx.oeEndB);
    EndForceResults{5,ii} = Chain4line.TimeHistory('End Force',ofx.Period(0, EndStageLength), ofx.oeEndB);
    EndForceResults{6,ii} = Chain4line.TimeHistory('End Force',ofx.Period(0, EndStageLength), ofx.oeEndA);
%     EndForceResults{7,ii} = Chain6line.TimeHistory('End Force',ofx.Period(0, EndStageLength), ofx.oeEndB);

    VonMisesResultsTime{1,ii} = Chain4line.TimeHistory('Max von Mises stress',ofx.Period(0, EndStageLength), ofx.oeEndB);


    MomentResults{1,ii} = MPline.RangeGraph('Bend Moment',ofx.arSpecifiedArclengths(0, TotalFoundationLength));
    VonMisesResults{1,ii} = MPline.RangeGraph('Max von Mises stress',ofx.arSpecifiedArclengths(0, TotalFoundationLength));
    ShearResults{1,ii} = MPline.RangeGraph('Shear force',ofx.arSpecifiedArclengths(0, TotalFoundationLength));
    DisplacementTowerResults{1,ii} = MPline.TimeHistory('X',ofx.Period(0, EndStageLength), ofx.oeEndB);    
    WaveElevationResults{1,ii} = environment.TimeHistory('Elevation', ofx.Period(0, EndStageLength), ofx.oeEnvironment(0, 0, 0));
    EffectiveTensionResults{1,ii} = MPline.RangeGraph('Effective tension',ofx.arSpecifiedArclengths(0, TotalFoundationLength));
    WallTensionResults{1,ii} = MPline.RangeGraph('Wall tension',ofx.arSpecifiedArclengths(0, TotalFoundationLength));
    ExternalPressureResults{1,ii} = MPline.RangeGraph('External pressure',ofx.arSpecifiedArclengths(0, TotalFoundationLength));
%     WholeSystemFrequencyResults1{1,ii} = 1/period1; % Find first frequency 1/period
%     WholeSystemFrequencyResults2{1,ii} = 1/period2; % Find second frequency 1/period
%     WholeSystemFrequencyResults3{1,ii} = 1/period3; % Find third frequency 1/period
    
    % From here cell values a extracted into arrays
    uMP(ii,:) = cell2mat(DisplacementResults(1,ii)); % 1 for MP, time history
    umaxMP(1,ii) = max(abs(cell2mat(DisplacementResults(1,ii)))); % 1 for MP, single value
    
    EndForce(ii,:) = cell2mat(EndForceResults(5,ii)); % 5 for chain 4 , time history
    EndForceMaxChain1(1,ii) = max(abs(cell2mat(EndForceResults(2,ii)))); % 2 for chain 1, single value
    EndForceMaxChain4(1,ii) = max(abs(cell2mat(EndForceResults(5,ii)))); % 5 for chain 4, single value
    EndForceMaxChain4EndA(1,ii) = max(abs(cell2mat(EndForceResults(6,ii)))); % 6 for chain 4 end A, single value
    EndForceMP(1,ii) = max((cell2mat(EndForceResults(1,ii)))); % 1 for MP, single value
    
    MomentMP(ii,:) = cell2mat(MomentResults(1,ii));
    VonMisesMP(ii,:) = cell2mat(VonMisesResults(1,ii));
    ShearMP(ii,:) = cell2mat(ShearResults(1,ii));
    uTower(ii,:) = cell2mat(DisplacementTowerResults(1,ii));
    umaxTower(1,ii) = max(abs(cell2mat(DisplacementTowerResults(1,ii))));
    eta(ii,:) = cell2mat(WaveElevationResults(1,ii));
    EffectiveTensionMP(ii,:) = cell2mat(EffectiveTensionResults(1,ii));
    WallTensionMP(ii,:) = cell2mat(WallTensionResults(1,ii));
    ExternalPressureMP(ii,:) = cell2mat(ExternalPressureResults(1,ii));
    
    VonMisesTime(ii,:) = cell2mat(VonMisesResultsTime(1,ii));
     
%     FreqSys1(ii,:) = cell2mat(WholeSystemFrequencyResults1(1,ii));
%     FreqSys2(ii,:) = cell2mat(WholeSystemFrequencyResults2(1,ii));
%     FreqSys3(ii,:) = cell2mat(WholeSystemFrequencyResults3(1,ii));

    l(1,ii) = LinesBC.EndAX(2); % Transfer footprint l values
    L(1,ii) = LinesBC.Length(2); % Transfer chain lengths L 
end % End geometry loop
    % Make sure all data from given input file are saved before being
    % overwritten 
    
    lvals(jj,:) = l;
    Lvals(jj,:) = L;
    umaxMPAllGeometry(jj,:) = umaxMP % All max displacements for all sheets for each input file jj    
    EndForceTimeHistory(jj,:) = EndForce;
    EndForceMaxAllGeometryChain1(jj,:) = EndForceMaxChain1 % All max end force for all sheets for each input file jj chain 1
    EndForceMaxAllGeometryChain4(jj,:) = EndForceMaxChain4 % All max end force for all sheets for each input file jj chain 1
    EndForceMaxAllGeometryChain4EndA(jj,:) = EndForceMaxChain4EndA % All max end force for all sheets for each input file jj chain 1
    EndForceMaxAllGeometryMP(jj,:) = EndForceMP % All max end force for all sheets for each input file jj chain 1
    MomentMPAllGeometry(jj,:) = MomentMP;
    VonMisesMPAllGeometry(jj,:) = VonMisesMP;
    ShearMPAllGeometry(jj,:) = ShearMP;
    umaxTowerAllGeometry(jj,:) = umaxTower;
    uTowerAllGeometry(jj,:) = uTower;
    etaAllGeometry(jj,:) = eta;
    EffectiveTensionMPAllGeometry(jj,:) = EffectiveTensionMP;
    WallTensionMPAllGeometry(jj,:) = WallTensionMP;
    ExternalPressureMPAllGeometry(jj,:) = ExternalPressureMP;
    VonMisesTimeAllGeomtery(jj,:) = VonMisesTime;
%     FreqSysAllGeometry1(jj,:) = FreqSys1
%     FreqSysAllGeometry2(jj,:) = FreqSys2
%     FreqSysAllGeometry3(jj,:) = FreqSys3
    fprintf('Geometry file finished \n')   
end % End NumberOfGeomInputFiles

end % End function


