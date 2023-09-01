function [AllOrbitTotalResult,AllOrbitAverageResult,ObjectMissionFlux,SigmResult] = ECP(Orbitdata,SimStartDates,Mode,ObjSpec,OrbitSimSpec)
    
    clc

    % Filepath for mastersettings
    ref_sdf = './Master/input/default.sdf';
    datei_sdf = '../Master/input/default.sdf';
    refpfad = './Master/input/master.inp';
    dateipfad = '../Master/input/master.inp';

    % Calculating sizes of all relevant surfaces
    [Surface,m] = DetermineSurfaceArea(ObjSpec); 

    % Initializing variables
    AllOrbitAverageResult = zeros(6,19,size(Orbitdata,1));
    AllOrbitTotalResult = zeros(6,19,size(Orbitdata,1));
    ObjectMissionFlux.TimeSections(1) = 0;
    SigmResult = zeros(6,4);
    
    % Clone Orbitdata and SimStartDates if more than one size and no atmosphere is set
    if strcmp(OrbitSimSpec.atmosphere,'None') && strcmp(ObjSpec.SpecOption,'Cubesat with Dragsail')
        Orbitdata{m,1} = NaN;
        for i=2 : m
            Orbitdata{i,1} = Orbitdata{1,1};
            SimStartDates(i) = SimStartDates(1);
        end
    end

    % Loop for size steps
    for i=1 : size(Orbitdata,1)     
        j=1;
        AktOrbit = Orbitdata{i,1};
        aktuellesDatum = SimStartDates(i);
        Result = zeros(6,19);
        Result(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
        Section = 0;
        
        % Mode 1 for complex mode with Master  
        if Mode == 1    

            % Set first heightstep
            SMANext = cell2mat(AktOrbit(1,6))-10;
            DeltatNext = cell2mat(AktOrbit(1,1))+(1/24);
            SimOrbits(1,:) = AktOrbit(1,:);
            k=2;

            % Evaluate all heightsteps
            while k <= size(AktOrbit, 1)
                if cell2mat(AktOrbit(k,1)) >= DeltatNext && cell2mat(AktOrbit(k,6)) < SMANext             
                    DeltatNext = DeltatNext+(1/24);
                    SMANext = SMANext-10;
                    j=j+1;
                    SimOrbits(j,:) = AktOrbit(k,:);
                end
                k=k+1;
            end

            % Loop for heightsteps
            for j=1 : size(SimOrbits, 1)     
                
                % Set identifier name
                Run_identifier = sprintf('Sphere                          -(27 char)-');
                
                % Set simulation date
                nDays = cell2mat(SimOrbits(j,1));
                SimBeginDate = aktuellesDatum + days(nDays);
                Year = year(SimBeginDate);
                Year = sprintf('%04d', Year);
                Month = month(SimBeginDate);
                Month = sprintf('%02d', Month);
                Day = day(SimBeginDate);
                Day = sprintf('%02d', Day);
                Hour = hour(SimBeginDate);
                Hour = sprintf('%02d', Hour);
                Sim_date_Begin = sprintf('%s %s %s %s  -(yyyy mm dd hh)- Begin',...
                    Year,...
                    Month,...
                    Day,...
                    Hour);
                Sim_date_End = sprintf('%s %s %s %s  -(yyyy mm dd hh)- End',...
                    Year,...
                    Month,...
                    Day,...
                    Hour);
                target_orbit = sprintf('0 %s %s %s %s %s %s %s %s %f %f %f %f %f', ...
                    Year, ...
                    Month, ...
                    Day, ...
                    Hour, ...
                    Year, ...
                    Month, ...
                    Day, ...
                    Hour, ...
                    string(SimOrbits(j,6)), ... % SMA
                    string(SimOrbits(j,7)), ... % ECC
                    string(SimOrbits(j,11)), ...% INC
                    string(SimOrbits(j,10)), ...% RAAN
                    string(SimOrbits(j,9)));    % AoP
                
                % Set comment
                SMA = cell2mat(SimOrbits(j,6));
                Run_comment = sprintf('Simulation with actual SMA %f',SMA);
                
                % Set target type
                if strcmp(ObjSpec.SpecOption,'Sphere')
                    target_type = '1               -(1:3)-            1 = sphere';    %1=Sphere 2=tumbling plate 3=oriented plate
                else
                    target_type = '3               -(1:3)-            1 = sphere';    %1=Sphere 2=tumbling plate 3=oriented plate
                end

                % Read input file
                fid = fopen(refpfad, 'r');
                inhalt = textscan(fid, '%s', 'delimiter', '\n');
                fclose(fid);
                
                % Change settings
                inhalt{1}{23} = Run_identifier;
                inhalt{1}{27} = Run_comment;
                inhalt{1}{34} = Sim_date_Begin;
                inhalt{1}{35} = Sim_date_End;
                inhalt{1}{99} = target_type;
                inhalt{1}{121} = target_orbit;
                
                % Write changed inputfile
                fid = fopen(refpfad, 'w');
                fprintf(fid, '%s\n', inhalt{1}{:});
                fclose(fid);
                
                % Copy input file to Masterpath
                copyfile(refpfad, dateipfad);
                copyfile(ref_sdf, datei_sdf); 

                % Starting Master
                
                % Specify the path to the executable
                executable_path = '..\Master\master-win64.exe';
                
                % Start the executable with the default application associated with its file type
                winopen(executable_path);
    
                % Check if application is still running
                while true
                    [~, result] = system('tasklist /fi "imagename eq master-win64.exe" /fo csv');
                    if contains(result, "master-win64.exe")
                       pause(1); 
                    else
                        break;
                    end
                end

                % Read solution file
                [AktData,SigmData] = MassData;
    
                % Evaluate time in heightstep
                if j == size(SimOrbits, 1)

                    % Set missionendtime
                    lastRow = AktOrbit{end, :};
    
                    % Calculate length of stay
                    nDays = lastRow - nDays;
                else

                    % Calculate length of stay
                    nDays = cell2mat(SimOrbits(j+1,1)) - nDays;
                end

                % Convert days to years
                t = nDays/365;
                
                % Collecting data for flux mission curve
                Section = Section + 1;
                SurfFlux.Sphere{i} = Surface.Sphere{1} .* (AktData./4);
                ObjectMissionFlux = ObjectMissionFluxCalc(SurfFlux,i,Section,nDays,ObjSpec.SpecOption,ObjectMissionFlux);

                % Multiply with length of stay
                AktData(:, 2:end) = AktData(:, 2:end) .* t;
                SigmData(:, 2:end) = SigmData(:, 2:end) .* t;

                % Add flux data 
                Result = Result + AktData; 

                % Add data for errorbars
                SigmResult = SigmResult + SigmData; 
            end

            % Collect results for all sizesteps
            AllOrbitAverageResult(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
            AllOrbitAverageResult(:, 2:end,1) = Result(:, 2:end) ./ (4*lastRow/365);
            SigmResult(:, 2:end) = SigmResult(:, 2:end) ./ (4*lastRow/365);
        
        % Mode 2 for reduced mode
        else
            nDaysBegin = 0;
            AktOrbit = cell2mat(AktOrbit);

            % Checking interval in line 1 and column 6
            SMAInterval = determine_interval(AktOrbit(1, 6), 10);

            % Checking interval in line 2 and column 11
            INCInterval = determine_interval(AktOrbit(1, 11), 10);  

            % Check timeinterval for first iteration

            % Evaluate remaining days of month
            maxAnzahlTage = eomday(year(aktuellesDatum), month(aktuellesDatum)); 
            vergangeneTage = getVergangeneTage(aktuellesDatum);
            DurationHr = maxAnzahlTage - vergangeneTage;  
            dauerInTagen = days(DurationHr);
            TimeInterval = double(dauerInTagen);
            
            DatabaseAltitudeInc = getDataBase(INCInterval,ObjSpec.SpecOption);
            SurfResult = CreateSurfResultStruct(ObjSpec.SpecOption,i);
            for k = 1:size(AktOrbit, 1)
                if AktOrbit(k, 6) <= SMAInterval-5 || AktOrbit(k,1) >= TimeInterval || k == size(AktOrbit, 1)

                    % Determine time interval
                    Section = Section + 1;
                    nDaysEnd = AktOrbit(k, 1);
                    nDaysCell = nDaysEnd - nDaysBegin;
                    nDaysBegin = nDaysEnd;                  
                    t = nDaysCell/365;

                    % Scale cell data                     
                    AktData = ScaleIncFlux(AktOrbit,INCInterval,SMAInterval,DatabaseAltitudeInc,ObjSpec.SpecOption);
                    AktData = ScaleTimeFlux(AktData,aktuellesDatum,AktOrbit(k,6),ObjSpec.SpecOption,AktOrbit(k,1),TimeInterval);
                    SurfFlux = ScaleSurfaceFlux(AktData,Surface,ObjSpec.SpecOption,i,AktOrbit(k, 6),AktOrbit(k, 11));

                    % Collecting data for flux mission curve
                    ObjectMissionFlux = ObjectMissionFluxCalc(SurfFlux,i,Section,nDaysEnd,ObjSpec.SpecOption,ObjectMissionFlux);

                    % Multiply with time of stay in cell
                    SurfFlux = ScaleTimeInterval(SurfFlux,t,ObjSpec.SpecOption,i);

                    % Add acutal value to result
                    SurfResult = AdaptResultSurfFlux(SurfResult,SurfFlux,ObjSpec.SpecOption,i);
                    
                    % check time interval
                    NextDate = aktuellesDatum + days(TimeInterval);
                    maxAnzahlTage = eomday(year(NextDate), month(NextDate));
                    TimeInterval = TimeInterval + maxAnzahlTage;

                    SMAInterval = determine_interval(AktOrbit(k, 6), 10);  % check SMA interval
                end
            end

            % Collect data for all sizesteps
            AllOrbitTotalResult = CollectAllOrbitTotalResult(SurfResult,ObjSpec.SpecOption,i,AllOrbitTotalResult);
            SurfResult = UnitCorrection(SurfResult,ObjSpec.SpecOption,nDaysEnd,i);
            AllOrbitAverageResult = CollectAllOrbitAverageResult(SurfResult,i,ObjSpec.SpecOption,Surface,AllOrbitAverageResult);
        end

    end
    if Mode == 1
        % implement save data function for complex mode
    else
        SaveData(AllOrbitAverageResult,AllOrbitTotalResult,ObjectMissionFlux,ObjSpec,OrbitSimSpec)
    end
end

function AllOrbitTotalResult = CollectAllOrbitTotalResult(SurfResult,SpecOption,i,AllOrbitTotalResult) % Collects the results of every orbit as total flux
    
    if strcmp(SpecOption,'Cubesat')
        % Assignment of fluxdata for the CubeSats
        TotalImpacts = zeros(6,19);
        TotalImpacts(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
        TotalImpacts(:, 2:end) = SurfResult.Cubesat.Front{i}(:, 2:end) + SurfResult.Cubesat.Right{i}(:, 2:end) + SurfResult.Cubesat.Left{i}(:, 2:end) + SurfResult.Cubesat.Up{i}(:, 2:end) + SurfResult.Cubesat.Down{i}(:, 2:end) + SurfResult.Cubesat.Rear{i}(:, 2:end);   
        AllOrbitTotalResult(:,:,i) = TotalImpacts;
    elseif strcmp(SpecOption,'Cubesat with Dragsail')
        % Assignment of fluxdata for the CubeSats and DragSailsteps
        TotalImpacts = zeros(6,19);
        TotalImpacts(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
        TotalImpacts(:, 2:end) = SurfResult.Cubesat.Front{i}(:, 2:end) + SurfResult.Cubesat.Right{i}(:, 2:end) + SurfResult.Cubesat.Left{i}(:, 2:end) + SurfResult.Cubesat.Up{i}(:, 2:end) + SurfResult.Cubesat.Down{i}(:, 2:end) + SurfResult.Cubesat.Rear{i}(:, 2:end); 
        AllOrbitTotalResult(:,:,i) = TotalImpacts;
        if i > 1
            TotalImpacts = zeros(6,19);
            TotalImpacts(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
            TotalImpacts(:, 2:end) = SurfResult.Dragsail.Front{i-1}(:, 2:end) + SurfResult.Dragsail.Rear{i-1}(:, 2:end)+SurfResult.Cubesat.Front{i}(:, 2:end) + SurfResult.Cubesat.Right{i}(:, 2:end) + SurfResult.Cubesat.Left{i}(:, 2:end) + SurfResult.Cubesat.Up{i}(:, 2:end) + SurfResult.Cubesat.Down{i}(:, 2:end) + SurfResult.Cubesat.Rear{i}(:, 2:end);   
            AllOrbitTotalResult(:,:,i) = TotalImpacts;
         
        end
    else
        AllOrbitTotalResult(:,:,i) = SurfResult.Sphere{i};
    end 
end

function SaveData(AllOrbitAverageResult,AllOrbitTotalResult,ObjectMissionFlux,ObjSpec,OrbitSimSpec) % Saves all data in a *.xlsx file 
    % Creating filename with date and time
    now = datetime('now','Format','dd-MMM-yyyy-HH-mm-ss');
    filename = fullfile('.','Output',sprintf('DebrisSim-%s.xlsx',now));

    sheetName = 'InputData'; % Name of Excel Workbook (Layer 3)
    if strcmp(ObjSpec.SpecOption,'Cubesat') 
        if strcmp(OrbitSimSpec.RefSys,'Orbital Elements') 
            PropData = {'Object Type', ObjSpec.SpecOption; ...
            'Width[cm]', ObjSpec.ObjWidth; ...
            'Height[cm]', ObjSpec.ObjHeight; ...
            'Depth[cm]', ObjSpec.ObjDepth; ...
            'Mass[kg]', ObjSpec.Mass; ...
            'Cd[]', ObjSpec.Cd; ...
            'Lift to Drag Ratio[]', ObjSpec.LtD;...
            '', ''; ...
            'Orbit Simulation Data', ''; ...
            'Orbital elements','';
            'SMA', OrbitSimSpec.a; ...
            'Ecc', OrbitSimSpec.e; ...
            'Inc', OrbitSimSpec.i; ...
            'LAN', OrbitSimSpec.LAN;...
            'AP', OrbitSimSpec.AP; ...
            'True Anomaly', OrbitSimSpec.f; ...
            'TL', OrbitSimSpec.TL; ...
            'AL', OrbitSimSpec.AL;...
            'LP', OrbitSimSpec.LP; ...
            '', ''; ...
            'Atmospheric model and Conditions','';...
            'Atmosphere Model', OrbitSimSpec.atmosphere;...
            'APH;', OrbitSimSpec.aph;...
            'Solar Flux', OrbitSimSpec.solarflux;...
            '', '';...
            'Gravity Model', OrbitSimSpec.GravityModel; ...
            '', '';...
            'Simulation Arguments','';...
            'ODE', OrbitSimSpec.ODE; ...
            'Simulation Time', OrbitSimSpec.td; ...
            'Final Height', OrbitSimSpec.hfinal;
            'Step Size', OrbitSimSpec.stepsize; ...
            'Rel. Error Tol.', OrbitSimSpec.RelTol; ...
            'Abs. Error Tol.', OrbitSimSpec.AbsTol;...
            '', '';...
            'Starttime', OrbitSimSpec.DateTime};
        else
             PropData = {'Object Type', ObjSpec.SpecOption; ...
            'Width[cm]', ObjSpec.ObjWidth; ...
            'Height[cm]', ObjSpec.ObjHeight; ...
            'Depth[cm]', ObjSpec.ObjDepth; ...
            'Mass[kg]', ObjSpec.Mass; ...
            'Cd[]', ObjSpec.Cd; ...
            'Lift to Drag Ratio[]', ObjSpec.LtD;
            '', ''; ...
            'Orbit Simulation Data', ''; ...
            'Spherical Coordinates','';
            'r0', OrbitSimSpec.r0; ...
            'h0', OrbitSimSpec.h0; ...
            'Phi0', OrbitSimSpec.Phi0; ...
            'Theta0', OrbitSimSpec.Theta0;...
            'rDot', OrbitSimSpec.rDot; ...
            'PhiDot', OrbitSimSpec.PhiDot; ...
            'ThetaDot', OrbitSimSpec.ThetaDot; ...
            '', ''; ...
            'Atmospheric model and Conditions','';
            'Atmosphere Model', OrbitSimSpec.atmosphere;...
            'APH;', OrbitSimSpec.aph;
            'Solar Flux', OrbitSimSpec.solarflux;
            '', '';
            'Gravity Model', OrbitSimSpec.GravityModel; ...
            '', '';
            'Simulation Arguments','';
            'ODE', OrbitSimSpec.ODE; ...
            'Simulation Time', OrbitSimSpec.td; ...
            'Final Height', OrbitSimSpec.hfinal;
            'Step Size', OrbitSimSpec.stepsize; ...
            'Rel. Error Tol.', OrbitSimSpec.RelTol; ...
            'Abs. Error Tol.', OrbitSimSpec.AbsTol;
            '', '';
            'Starttime', OrbitSimSpec.DateTime};
        end
    elseif strcmp(ObjSpec.SpecOption,'Cubesat with Dragsail')
        if strcmp(OrbitSimSpec.RefSys,'Orbital Elements')
            PropData = {'Object Type', ObjSpec.SpecOption; ...
            'Width[cm]', ObjSpec.ObjWidth; ...
            'Height[cm]', ObjSpec.ObjHeight; ...
            'Depth[cm]', ObjSpec.ObjDepth; ...
            'MinProjDragsailArea[m^2]', ObjSpec.MinProjDragsailArea; ...
            'MaxProjDragsailArea[m^2]', ObjSpec.MaxProjDragsailArea; ...
            'Stepnumber[]', ObjSpec.Stepnumber; ...
            'Mass[kg]', ObjSpec.Mass; ...
            'Cd[]', ObjSpec.Cd; ...
            'Lift to Drag Ratio[]', ObjSpec.LtD; 
            '', ''; ...
            'Orbit Simulation Data', ''; ...
            'Orbital elements','';
            'SMA', OrbitSimSpec.a; ...
            'Ecc', OrbitSimSpec.e; ...
            'Inc', OrbitSimSpec.i; ...
            'LAN', OrbitSimSpec.LAN;...
            'AP', OrbitSimSpec.AP; ...
            'True Anomaly', OrbitSimSpec.f; ...
            'TL', OrbitSimSpec.TL; ...
            'AL', OrbitSimSpec.AL;
            'LP', OrbitSimSpec.LP; ...
            '', ''; ...
            'Atmospheric model and Conditions','';
            'Atmosphere Model', OrbitSimSpec.atmosphere;...
            'APH;', OrbitSimSpec.aph;
            'Solar Flux', OrbitSimSpec.solarflux;
            '', '';
            'Gravity Model', OrbitSimSpec.GravityModel; ...
            '', '';
            'Simulation Arguments','';
            'ODE', OrbitSimSpec.ODE; ...
            'Simulation Time', OrbitSimSpec.td; ...
            'Final Height', OrbitSimSpec.hfinal;
            'Step Size', OrbitSimSpec.stepsize; ...
            'Rel. Error Tol.', OrbitSimSpec.RelTol; ...
            'Abs. Error Tol.', OrbitSimSpec.AbsTol;
            '', '';
            'Starttime', OrbitSimSpec.DateTime};
        else
            PropData = {'Object Type', ObjSpec.SpecOption; ...
            'Width[cm]', ObjSpec.ObjWidth; ...
            'Height[cm]', ObjSpec.ObjHeight; ...
            'Depth[cm]', ObjSpec.ObjDepth; ...
            'MinProjDragsailArea[m^2]', ObjSpec.MinProjDragsailArea; ...
            'MaxProjDragsailArea[m^2]', ObjSpec.MaxProjDragsailArea; ...
            'Stepnumber[]', ObjSpec.Stepnumber; ...
            'Mass[kg]', ObjSpec.Mass; ...
            'Cd[]', ObjSpec.Cd; ...
            'Lift to Drag Ratio[]', ObjSpec.LtD; 
            '', ''; ...
            'Orbit Simulation Data', ''; ...
            'Spherical Coordinates','';
            'r0', OrbitSimSpec.r0; ...
            'h0', OrbitSimSpec.h0; ...
            'Phi0', OrbitSimSpec.Phi0; ...
            'Theta0', OrbitSimSpec.Theta0;...
            'rDot', OrbitSimSpec.rDot; ...
            'PhiDot', OrbitSimSpec.PhiDot; ...
            'ThetaDot', OrbitSimSpec.ThetaDot; ...
            '', ''; ...
            'Atmospheric model and Conditions','';
            'Atmosphere Model', OrbitSimSpec.atmosphere;...
            'APH;', OrbitSimSpec.aph;
            'Solar Flux', OrbitSimSpec.solarflux;
            '', '';
            'Gravity Model', OrbitSimSpec.GravityModel; ...
            '', '';
            'Simulation Arguments','';
            'ODE', OrbitSimSpec.ODE; ...
            'Simulation Time', OrbitSimSpec.td; ...
            'Final Height', OrbitSimSpec.hfinal;
            'Step Size', OrbitSimSpec.stepsize; ...
            'Rel. Error Tol.', OrbitSimSpec.RelTol; ...
            'Abs. Error Tol.', OrbitSimSpec.AbsTol;
            '', '';
            'Starttime', OrbitSimSpec.DateTime};
        end
    else
        if strcmp(OrbitSimSpec.RefSys,'Orbital Elements')
            PropData = {'Object Data', '';...
            'Object Type', ObjSpec.SpecOption; ...
            'ballistic coeficient[kg/m^2]', ObjSpec.bInput; ...
            'Mass[kg]', ObjSpec.Mass; ...
            'Cd[]', ObjSpec.Cd; ...
            'Lift to Drag ratio[]', ObjSpec.LtD;...
            '', ''; ...
            'Orbit Simulation Data', ''; ...
            'Orbital elements','';
            'SMA', OrbitSimSpec.a; ...
            'Ecc', OrbitSimSpec.e; ...
            'Inc', OrbitSimSpec.i; ...
            'LAN', OrbitSimSpec.LAN;...
            'AP', OrbitSimSpec.AP; ...
            'True Anomaly', OrbitSimSpec.f; ...
            'TL', OrbitSimSpec.TL; ...
            'AL', OrbitSimSpec.AL;
            'LP', OrbitSimSpec.LP; ...
            '', ''; ...
            'Atmospheric model and Conditions','';
            'Atmosphere Model', OrbitSimSpec.atmosphere;...
            'APH;', OrbitSimSpec.aph;
            'Solar Flux', OrbitSimSpec.solarflux;
            '', '';
            'Gravity Model', OrbitSimSpec.GravityModel; ...
            '', '';
            'Simulation Arguments','';
            'ODE', OrbitSimSpec.ODE; ...
            'Simulation Time', OrbitSimSpec.td; ...
            'Final Height', OrbitSimSpec.hfinal;
            'Step Size', OrbitSimSpec.stepsize; ...
            'Rel. Error Tol.', OrbitSimSpec.RelTol; ...
            'Abs. Error Tol.', OrbitSimSpec.AbsTol;
            '', '';
            'Starttime', OrbitSimSpec.DateTime};
        else
            PropData = {'Object Data', '';...
            'Object Type', ObjSpec.SpecOption; ...
            'ballistic coeficient[kg/m^2]', ObjSpec.bInput; ...
            'Mass[kg]', ObjSpec.Mass; ...
            'Cd[]', ObjSpec.Cd; ...
            'Lift to Drag ratio[]', ObjSpec.LtD;...
            '', ''; ...
            'Orbit Simulation Data', ''; ...
            'Spherical Coordinates','';
            'r0', OrbitSimSpec.r0; ...
            'h0', OrbitSimSpec.h0; ...
            'Phi0', OrbitSimSpec.Phi0; ...
            'Theta0', OrbitSimSpec.Theta0;...
            'rDot', OrbitSimSpec.rDot; ...
            'PhiDot', OrbitSimSpec.PhiDot; ...
            'ThetaDot', OrbitSimSpec.ThetaDot; ...
            '', ''; ...
            'Atmospheric model and Conditions','';
            'Atmosphere Model', OrbitSimSpec.atmosphere;...
            'APH;', OrbitSimSpec.aph;
            'Solar Flux', OrbitSimSpec.solarflux;
            '', '';
            'Gravity Model', OrbitSimSpec.GravityModel; ...
            '', '';
            'Simulation Arguments','';
            'ODE', OrbitSimSpec.ODE; ...
            'Simulation Time', OrbitSimSpec.td; ...
            'Final Height', OrbitSimSpec.hfinal;
            'Step Size', OrbitSimSpec.stepsize; ...
            'Rel. Error Tol.', OrbitSimSpec.RelTol; ...
            'Abs. Error Tol.', OrbitSimSpec.AbsTol;
            '', '';
            'Starttime', OrbitSimSpec.DateTime};
        end
    end

    writecell(PropData, filename, 'Sheet', sheetName);
    % Open excelfile with COM-objects
    excelObj = actxserver('Excel.Application');
    workbookObj = excelObj.Workbooks.Open(fullfile(pwd, filename));
    sheetObj = workbookObj.Sheets.Item(sheetName);
    
    % Change range of the second column
    columnBRange = sheetObj.Range('B:B');
    columnBRange.HorizontalAlignment = -4131;  % Konstant for align left
    
    % adapt column width
    columnBRange.EntireColumn.AutoFit();

    % Save Excelfile and close
    workbookObj.Save();
    workbookObj.Close();
    excelObj.Quit();

    sheetName = 'ObjectMissionFlux per year'; % Name of workbook 
    if strcmp(ObjSpec.SpecOption,'Cubesat') 
        % Ihre Originaldatenmatrix
        originalData = ObjectMissionFlux.Cubesat.Total;
        
        % Die Größe der Originaldatenmatrix
        [rows, ~] = size(originalData);
        
        % Entfernen jeder 8. Zeile
        indexToRemove = 8:8:rows;
        originalData(indexToRemove, :) = [];
        
        % Neue Größe der Datenmatrix
        [rows, cols] = size(originalData);
        
        % Die Anzahl der zusätzlichen Zeilen für die Strings
        numExtraRows = floor(rows / 7);  % +1 für den ersten Titel
        
        % Erstellen einer leeren Zellenmatrix
        CubesatData = cell(rows + numExtraRows, cols);
        
        % Einfügen der Originaldaten und der Strings
        count = 1;
        titleCount = 1;
        
        for i = 1:(rows + numExtraRows)
            if i == 1 || mod(i - 1, 8) == 0
                CubesatData{i, 1} = ['Cubesat' num2str(titleCount-1)];
                titleCount = titleCount + 1;
            else
                CubesatData(i, :) = num2cell(originalData(count, :));
                count = count + 1;
            end
        end

        data = CubesatData;
    elseif strcmp(ObjSpec.SpecOption,'Cubesat with Dragsail')
        % Ihre Originaldatenmatrix
        originalData = ObjectMissionFlux.Cubesat.Total;
        
        % Die Größe der Originaldatenmatrix
        [rows, ~] = size(originalData);
        
        % Entfernen jeder 8. Zeile
        indexToRemove = 8:8:rows;
        originalData(indexToRemove, :) = [];
        
        % Neue Größe der Datenmatrix
        [rows, cols] = size(originalData);
        
        % Die Anzahl der zusätzlichen Zeilen für die Strings
        numExtraRows = floor(rows / 7);  % +1 für den ersten Titel
        
        % Erstellen einer leeren Zellenmatrix
        CubesatData = cell(rows + numExtraRows, cols);
        
        % Einfügen der Originaldaten und der Strings
        count = 1;
        titleCount = 1;
        
        for i = 1:(rows + numExtraRows)
            if i == 1 || mod(i - 1, 8) == 0
                CubesatData{i, 1} = ['Cubesat' num2str(titleCount-1)];
                titleCount = titleCount + 1;
            else
                CubesatData(i, :) = num2cell(originalData(count, :));
                count = count + 1;
            end
        end   
    
        % Ihre Originaldatenmatrix
        originalData = ObjectMissionFlux.Dragsail.Total;
        
        % Die Größe der Originaldatenmatrix
        [rows, ~] = size(originalData);
        
        % Entfernen jeder 8. Zeile
        indexToRemove = 8:8:rows;
        originalData(indexToRemove, :) = [];
        
        % Neue Größe der Datenmatrix
        [rows, cols] = size(originalData);
        
        % Die Anzahl der zusätzlichen Zeilen für die Strings
        numExtraRows = floor(rows / 7);  % +1 für den ersten Titel
        
        % Erstellen einer leeren Zellenmatrix
        DragsailData = cell(rows + numExtraRows, cols);
        
        % Einfügen der Originaldaten und der Strings
        count = 1;
        titleCount = 1;
        
        for i = 1:(rows + numExtraRows)
            if i == 1 || mod(i - 1, 8) == 0
                DragsailData{i, 1} = ['DragSail' num2str(titleCount)];
                titleCount = titleCount + 1;
            else
                DragsailData(i, :) = num2cell(originalData(count, :));
                count = count + 1;
            end
        end     

        data = [CubesatData ; DragsailData];
    else
        % Ihre Originaldatenmatrix
        originalData = ObjectMissionFlux.Sphere;
        
        % Die Größe der Originaldatenmatrix
        [rows, ~] = size(originalData);
        
        % Entfernen jeder 8. Zeile
        indexToRemove = 8:8:rows;
        originalData(indexToRemove, :) = [];
        
        % Neue Größe der Datenmatrix
        [rows, cols] = size(originalData);
        
        % Die Anzahl der zusätzlichen Zeilen für die Strings
        numExtraRows = floor(rows / 7);  % +1 für den ersten Titel
        
        % Erstellen einer leeren Zellenmatrix
        SphereData = cell(rows + numExtraRows, cols);
        
        % Einfügen der Originaldaten und der Strings
        count = 1;
        titleCount = 1;
        
        for i = 1:(rows + numExtraRows)
            if i == 1 || mod(i - 1, 8) == 0
                SphereData{i, 1} = ['Sphere' num2str(titleCount-1)];
                titleCount = titleCount + 1;
            else
                SphereData(i, :) = num2cell(originalData(count, :));
                count = count + 1;
            end
        end

        data = SphereData;
    end
    
    % Zeilentitel erstellen
    rowTitles = {'','Time [d]','10^-3 bis 10^-2 [kg]','10^-2 bis 10^-1 [kg]','10^-1 bis 10^0 [kg]','10^0 bis 10^1 [kg]','10^1 bis 10^2 [kg]','10^2 bis 10^3 [kg]'}; % Mit Leerzeile am Ende
    numRows = size(data, 1);
    numTitleRows = length(rowTitles);
    fullRowTitles = repmat(rowTitles, 1, ceil(numRows/numTitleRows));
    fullRowTitles = fullRowTitles(1:numRows); % Auf die Gesamtzeilenzahl der Matrix kürzen
    
    % Daten und Titel zusammenfügen
    cellData = cell(size(data, 1), size(data, 2) + 1);
    cellData(:,1) = fullRowTitles(:);
    cellData(:,2:end) = data;
    
    % Daten in Excel schreiben, beginnend bei der ersten Spalte (A)
    writecell(cellData, filename, 'Sheet', sheetName, 'Range', 'A1');

    % Write average results in file
    for i = 1:size(AllOrbitAverageResult, 3)
        if strcmp(ObjSpec.SpecOption,'Cubesat with Dragsail')
            if i == 1
                sheetName = ['ObjectAverageFlux_CubeSat', num2str(i-1)]; % Name of workbook 
            else
                sheetName = ['ObjectAverageFlux_Dragsail', num2str(i-1)]; % Name of workbook 
            end
        else
            sheetName = ['ObjectAverageFlux', num2str(i)]; % Name of workbook 
        end
        data = AllOrbitAverageResult(:, :, i);
        writematrix(data, filename, 'Sheet', sheetName);
    end
    
    % Write total results in file
    for i = 1:size(AllOrbitTotalResult, 3)
        if strcmp(ObjSpec.SpecOption,'Cubesat with Dragsail')
            if i == 1
                sheetName = ['ObjectTotalFlux_CubeSat', num2str(i-1)]; % Name of workbook 
            else
                sheetName = ['ObjectTotalFlux_Dragsail', num2str(i-1)]; % Name of workbook 
            end
        else
            sheetName = ['ObjectTotalFlux', num2str(i)]; % Name of workbook 
        end
        data = AllOrbitTotalResult(:, :, i);
        writematrix(data, filename, 'Sheet', sheetName);
    end
end

function ObjectMissionFlux = ObjectMissionFluxCalc(SurfFlux,i,Section,nDaysEnd,SpecOption,ObjectMissionFlux) % Collects the results of every CPE to provide flux history 
    ObjectMissionFlux.TimeSections(i,Section) = nDaysEnd;
    if strcmp(SpecOption,'Cubesat')
        % Assignment of fluxdata for CubeSats
        ObjectMissionFlux.Cubesat.Front(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Front(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Front{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Right(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Right(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Right{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Left(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Left(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Left{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Up(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Up(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Up{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Down(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Down(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Down{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Rear(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Rear(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Rear{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Total(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Total(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Rear{i}(1:6, 19) + SurfFlux.Cubesat.Down{i}(1:6, 19) + SurfFlux.Cubesat.Up{i}(1:6, 19) + SurfFlux.Cubesat.Left{i}(1:6, 19) + SurfFlux.Cubesat.Right{i}(1:6, 19) + SurfFlux.Cubesat.Front{i}(1:6, 19);
    elseif strcmp(SpecOption,'Cubesat with Dragsail')
        % Assignment of fluxdata for CubeSats and Dragsailsteps
        
        ObjectMissionFlux.Cubesat.Front(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Front(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Front{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Right(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Right(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Right{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Left(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Left(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Left{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Up(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Up(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Up{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Down(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Down(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Down{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Rear(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Rear(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Rear{i}(1:6, 19);
        ObjectMissionFlux.Cubesat.Total(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Cubesat.Total(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Cubesat.Rear{i}(1:6, 19) + SurfFlux.Cubesat.Down{i}(1:6, 19) + SurfFlux.Cubesat.Up{i}(1:6, 19) + SurfFlux.Cubesat.Left{i}(1:6, 19) + SurfFlux.Cubesat.Right{i}(1:6, 19) + SurfFlux.Cubesat.Front{i}(1:6, 19);
        if i > 1 
            ObjectMissionFlux.Dragsail.Front(1+(8*(i-2)),Section) = nDaysEnd;
            ObjectMissionFlux.Dragsail.Front(2+(8*(i-2)):7+(8*(i-2)),Section) = SurfFlux.Dragsail.Front{i-1}(1:6, 19) ;
            ObjectMissionFlux.Dragsail.Rear(1+(8*(i-2)),Section) = nDaysEnd;
            ObjectMissionFlux.Dragsail.Rear(2+(8*(i-2)):7+(8*(i-2)),Section) = SurfFlux.Dragsail.Rear{i-1}(1:6, 19) ;
            ObjectMissionFlux.Dragsail.Total(1+(8*(i-2)),Section) = nDaysEnd;
            ObjectMissionFlux.Dragsail.Total(2+(8*(i-2)):7+(8*(i-2)),Section) = SurfFlux.Dragsail.Front{i-1}(1:6, 19) + SurfFlux.Dragsail.Front{i-1}(1:6, 19)+SurfFlux.Cubesat.Rear{i}(1:6, 19) + SurfFlux.Cubesat.Down{i}(1:6, 19) + SurfFlux.Cubesat.Up{i}(1:6, 19) + SurfFlux.Cubesat.Left{i}(1:6, 19) + SurfFlux.Cubesat.Right{i}(1:6, 19) + SurfFlux.Cubesat.Front{i}(1:6, 19);
        end
    else
        ObjectMissionFlux.Sphere(1+(8*(i-1)),Section) = nDaysEnd;
        ObjectMissionFlux.Sphere(2+(8*(i-1)):7+(8*(i-1)),Section) = SurfFlux.Sphere{i}(1:6, 19);
    end
end

function AllOrbitAverageResult = CollectAllOrbitAverageResult(SurfResult,i,SpecOption,Surface,AllOrbitAverageResult) % Collects the results of ervery orbit as average flux
    
    if strcmp(SpecOption,'Cubesat')
        % Assignment of fluxdata for CubeSats
        AverageFlux = zeros(6,19);
        AverageFlux(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
        AverageFlux(:, 2:end) = SurfResult.Cubesat.Front{i}(:, 2:end) + SurfResult.Cubesat.Right{i}(:, 2:end) + SurfResult.Cubesat.Left{i}(:, 2:end) + SurfResult.Cubesat.Up{i}(:, 2:end) + SurfResult.Cubesat.Down{i}(:, 2:end) + SurfResult.Cubesat.Rear{i}(:, 2:end);   
        AverageFlux(:, 2:end) =AverageFlux(:, 2:end) ./ (Surface.Cubesat.Front + Surface.Cubesat.Up + Surface.Cubesat.Down + Surface.Cubesat.Right + Surface.Cubesat.Left + Surface.Cubesat.Rear);
        AllOrbitAverageResult(:,:,i) = AverageFlux;
    elseif strcmp(SpecOption,'Cubesat with Dragsail')
        % Assignment of fluxdata for CubeSats and DragSailsteps

        AverageFlux = zeros(6,19);
        AverageFlux(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
        AverageFlux(:, 2:end) = SurfResult.Cubesat.Front{i}(:, 2:end) + SurfResult.Cubesat.Right{i}(:, 2:end) + SurfResult.Cubesat.Left{i}(:, 2:end) + SurfResult.Cubesat.Up{i}(:, 2:end) + SurfResult.Cubesat.Down{i}(:, 2:end) + SurfResult.Cubesat.Rear{i}(:, 2:end); 
        AverageFlux(:, 2:end) =AverageFlux(:, 2:end) ./ (Surface.Cubesat.Front + Surface.Cubesat.Up + Surface.Cubesat.Down + Surface.Cubesat.Right + Surface.Cubesat.Left + Surface.Cubesat.Rear);
        AllOrbitAverageResult(:,:,i) = AverageFlux;
        if i > 1 
                AverageFlux = zeros(6,19);
                AverageFlux(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
                AverageFlux(:, 2:end) = SurfResult.Dragsail.Front{i-1}(:, 2:end) + SurfResult.Dragsail.Rear{i-1}(:, 2:end)+SurfResult.Cubesat.Front{i}(:, 2:end) + SurfResult.Cubesat.Right{i}(:, 2:end) + SurfResult.Cubesat.Left{i}(:, 2:end) + SurfResult.Cubesat.Up{i}(:, 2:end) + SurfResult.Cubesat.Down{i}(:, 2:end) + SurfResult.Cubesat.Rear{i}(:, 2:end); 
                AverageFlux(:, 2:end) =AverageFlux(:, 2:end) ./ (Surface.Dragsail.Front{i-1}+Surface.Dragsail.Rear{i-1}+Surface.Cubesat.Front + Surface.Cubesat.Up + Surface.Cubesat.Down + Surface.Cubesat.Right + Surface.Cubesat.Left + Surface.Cubesat.Rear);     
                AllOrbitAverageResult(:,:,i) = AverageFlux;
        end
    else
        SurfResult.Sphere{i}(:, 2:end) = SurfResult.Sphere{i}(:, 2:end) ./Surface.Sphere{1};
        AllOrbitAverageResult(:,:,i) = SurfResult.Sphere{i};
    end 
end

function vergangeneTage = getVergangeneTage(aktuellesDatum) % extracts the passed days and hours of a month
    aktuelleStunde = hour(aktuellesDatum); % present hour
    vergangeneTage = day(aktuellesDatum) + hours(aktuelleStunde) -1; % present hour and day
end

function interval = determine_interval(value, interval_size) % Function for determining the Cellinterval
    interval_start = floor(value / interval_size) * interval_size;
    interval_end = interval_start + interval_size;
    interval_midpoint = (interval_start + interval_end) / 2;
    interval = interval_midpoint;
end

function Interval = GetInterval(Value,startwert,Intervalsize) % finds the interval for the database entries
    
    naechsterWert = startwert;
    differenz = abs(Value - naechsterWert);
    
    while true
        naechsterWert = naechsterWert + Intervalsize; % Calculates the next interval 
        neueDifferenz = abs(Value - naechsterWert);
        
        if neueDifferenz < differenz
            differenz = neueDifferenz;
        else
            Interval = naechsterWert - Intervalsize;
            break; % Stop loop, if difference is bigger or equal
        end
    end
    if Interval > 8375
        Interval = 8375;
    end
end

function  DatabaseAltitudeInc = getDataBase(INCInterval,SpecOption) % Gets the correct database file

    % Filename for Excelfile
    if strcmp(SpecOption,'Sphere')
        filename_excel = sprintf('Impact-Flux-Database/FluxAltitudeInclination/FluxAltitudeMass_INC%d.xlsx',INCInterval);
    else
        filename_excel = sprintf('Impact-Flux-Database/FluxAltitudeInclination_FS/FluxAltitudeMass_FS_INC%d.xlsx',INCInterval);
    end
    
    % Read Excelfile
    DatabaseAltitudeInc = readmatrix(filename_excel);
end

function  DatabaseTimeAltitude = getTimeDataBase(SMATimeInterval,SpecOption) % Gets the correct database file
    % Filename for Excelfile
    if strcmp(SpecOption,'Sphere')
        filename_excel = sprintf('Impact-Flux-Database/FluxTimeAltitude/FluxTimeMass_Altitude%d.xlsx',SMATimeInterval);
    else
        filename_excel = sprintf('Impact-Flux-Database/FluxTimeAltitude_FS/FluxTimeMass_FS_Altitude%d.xlsx',SMATimeInterval);
    end
    
    % Read Excelfile
    DatabaseTimeAltitude = readmatrix(filename_excel);
end

function [AktData,SigmData] = MassData() % Get MASTER data
    % Current rownumber
    zeilennummer = 1;

    % Filename of datafile
    dateiname = '..\Master\output\Sphere_d.mas';
    
    % Open file
    datei = fopen(dateiname, 'r');
    
    % Initialize table
    AktData = zeros(6,19);
    
    % Read row data
    zeile = fgetl(datei);
    while ischar(zeile)
        % Check if line starts with "#" 
        if ~startsWith(zeile, '#')
            % Split each line at spaces or tabs
            zeilendaten = strsplit(zeile, {' ', '\t'});
            
            % Remove empty elements
            zeilendaten = zeilendaten(~cellfun('isempty', zeilendaten));
            
            % Convert row data in numerical values
            zeilendaten = str2double(zeilendaten);
            
            % Check if row data is valid
            if ~any(isnan(zeilendaten))
                % Add row data to table
                AktData(zeilennummer, :) = zeilendaten;
                
                % Inkrementiere die Zeilennummer
                zeilennummer = zeilennummer + 1;
            end
        end
        
        % Read next line
        zeile = fgetl(datei);
    end
    
    % Close file
    fclose(datei);

    %_______________________Uncertainty bars_____________________________
    
    % Current rownumber
    zeilennummer = 1;    

    % Filename of datafile
    dateiname = '..\Master\output\Sphere_d.mas.sigm';
    
    % Open file
    datei = fopen(dateiname, 'r');
    
    % Initialize table
    SigmData = zeros(6,4);
    
    % Read row data
    zeile = fgetl(datei);
    while ischar(zeile)
        % Check if line starts with "#"
        if ~startsWith(zeile, '#')
            % Split each line at spaces or tabs
            zeilendaten = strsplit(zeile, {' ', '\t'});
            
            % Remove empty elements
            zeilendaten = zeilendaten(~cellfun('isempty', zeilendaten));
            
            % Convert row data into numerical values
            zeilendaten = str2double(zeilendaten);
            
            % Check if row data is valid
            if ~any(isnan(zeilendaten))
                % Add row data to table
                SigmData(zeilennummer, :) = zeilendaten;
                
                % Increment the row number
                zeilennummer = zeilennummer + 1;
            end
        end
        
        % Read next line
        zeile = fgetl(datei);
    end
    
    % Close file
    fclose(datei);
end

function [Surface,m] = DetermineSurfaceArea(ObjSpec) % Calculates the surface area of the satellite with the input values
    if strcmp(ObjSpec.SpecOption,'Cubesat')
        % Calculate all cubeSats surfaces
        m = 1;
        Surface.Cubesat.Front = ObjSpec.ObjWidth*1e-2*ObjSpec.ObjHeight*1e-2;
        Surface.Cubesat.Right = ObjSpec.ObjDepth*1e-2*ObjSpec.ObjHeight*1e-2;
        Surface.Cubesat.Left = ObjSpec.ObjDepth*1e-2*ObjSpec.ObjHeight*1e-2;
        Surface.Cubesat.Up = ObjSpec.ObjWidth*1e-2*ObjSpec.ObjDepth*1e-2;
        Surface.Cubesat.Down = ObjSpec.ObjWidth*1e-2*ObjSpec.ObjDepth*1e-2;
        Surface.Cubesat.Rear = ObjSpec.ObjWidth*1e-2*ObjSpec.ObjHeight*1e-2;
    elseif strcmp(ObjSpec.SpecOption,'Cubesat with Dragsail')
        % Calculate all cubeSats and dragSail surfaces
        m = 1 + ObjSpec.Stepnumber;
        Surface.Cubesat.Front = ObjSpec.ObjWidth*1e-2*ObjSpec.ObjHeight*1e-2;
        Surface.Cubesat.Right = ObjSpec.ObjDepth*1e-2*ObjSpec.ObjHeight*1e-2;
        Surface.Cubesat.Left = ObjSpec.ObjDepth*1e-2*ObjSpec.ObjHeight*1e-2;
        Surface.Cubesat.Up = ObjSpec.ObjWidth*1e-2*ObjSpec.ObjDepth*1e-2;
        Surface.Cubesat.Down = ObjSpec.ObjWidth*1e-2*ObjSpec.ObjDepth*1e-2;
        Surface.Cubesat.Rear = ObjSpec.ObjWidth*1e-2*ObjSpec.ObjHeight*1e-2;
        if ObjSpec.Stepnumber == 1
            Surface.Dragsail.Front{1} = ObjSpec.MinProjDragsailArea;
            Surface.Dragsail.Rear{1} = ObjSpec.MinProjDragsailArea;
        elseif ObjSpec.Stepnumber == 2
            Surface.Dragsail.Front{1} = ObjSpec.MinProjDragsailArea;
            Surface.Dragsail.Rear{1} = ObjSpec.MinProjDragsailArea;
            Surface.Dragsail.Front{2} = ObjSpec.MaxProjDragsailArea;
            Surface.Dragsail.Rear{2} = ObjSpec.MaxProjDragsailArea;
        elseif ObjSpec.Stepnumber > 2
            Surface.Dragsail.Front{1} = ObjSpec.MinProjDragsailArea;
            Surface.Dragsail.Rear{1} = ObjSpec.MinProjDragsailArea;
            
            for j=3 : ObjSpec.Stepnumber
                Stepsize = (ObjSpec.MaxProjDragsailArea-ObjSpec.MinProjDragsailArea)/(ObjSpec.Stepnumber-1);
                Surface.Dragsail.Front{j-1} = ObjSpec.MinProjDragsailArea+(j-2)*Stepsize;
                Surface.Dragsail.Rear{j-1} = ObjSpec.MinProjDragsailArea+(j-2)*Stepsize;
            end
            Surface.Dragsail.Front{ObjSpec.Stepnumber} = ObjSpec.MaxProjDragsailArea;
            Surface.Dragsail.Rear{ObjSpec.Stepnumber} = ObjSpec.MaxProjDragsailArea;
        end
    else
        % Calculating crosssection out of ballistic coeficient    
        b = strsplit(ObjSpec.bInput); % Converts string of ballistic coefficient to cell array
        b =str2double(b); % Converts strings to numeric values
        for j=1 : size(b,2) % Loop for different ballistic coefficients 
        Surface.Sphere{j} = (ObjSpec.Mass / ObjSpec.Cd *b(j))*4;
        end
        m = size(b,2);
    end
end

function AktData = ScaleIncFlux(AktOrbit,INCInterval,SMAInterval,DatabaseAltitudeIncActCell,SpecOption) % Scales the flux data to the correct value for the current inclination
    if AktOrbit(end,11) > INCInterval 
        if AktOrbit(end,11) > 175
            [row, ~] = find(DatabaseAltitudeIncActCell == SMAInterval);
            AktData(:, 19)=DatabaseAltitudeIncActCell(row, 2:end);
        else
            INCIntervalNext = INCInterval + 10;
            DatabaseAltitudeInc = getDataBase(INCIntervalNext,SpecOption);
            [row, ~] = find(DatabaseAltitudeInc == SMAInterval);
            if isempty(row)
                row = size(DatabaseAltitudeInc,1);
            end
            AktData(:, 19)=(DatabaseAltitudeIncActCell(row, 2:end)+((DatabaseAltitudeInc(row, 2:end)-DatabaseAltitudeIncActCell(row, 2:end))*(AktOrbit(end,11)-INCInterval))/(INCIntervalNext-INCInterval));     
        end
    else
        if AktOrbit(end,11) <= 5
            [row, ~] = find(DatabaseAltitudeIncActCell == SMAInterval);
            AktData(:, 19)=DatabaseAltitudeIncActCell(row, 2:end);
        else
            INCIntervalNext = INCInterval - 10;
            DatabaseAltitudeInc = getDataBase(INCIntervalNext,SpecOption);
            [row, ~] = find(DatabaseAltitudeInc == SMAInterval);
            if isempty(row)
                row = size(DatabaseAltitudeInc,1);
            end
            AktData(:, 19)=(DatabaseAltitudeIncActCell(row, 2:end)-((DatabaseAltitudeInc(row, 2:end)-DatabaseAltitudeIncActCell(row, 2:end))*(AktOrbit(end,11)-INCInterval))/(INCIntervalNext-INCInterval));
        end
    end 

end

function AktData = ScaleTimeFlux(AktData,startDatum,AktSMA,SpecOption,AktTime,TimeIntervalTop) % Adapts the flux data for different epochs
    % Convert the start date to a MATLAB date format
    startDatum = datetime(startDatum);
    
    % Calculates the final date
    endDatum = startDatum + days(AktTime);
 
    % Create the reference date (1. January 2023)
    referenzDatum = datetime(2023, 1, 1);
    
    % Calculate the month difference
    TimeDuration = between(referenzDatum,endDatum,"months");
    TimeInterval = calmonths(TimeDuration); 
    if AktTime < TimeIntervalTop
        TimeInterval = TimeInterval+1;
    elseif TimeInterval < 1 % If the date is before the reference month, the function will exit and not scale
        return
    end

    % Determine the SMA interval to find the correct file
    startwert = 6375;
    Intervalsize = 100;
    SMATimeInterval = GetInterval(AktSMA,startwert,Intervalsize); % Finding interval
    DatabaseTimeAltitude = getTimeDataBase(SMATimeInterval,SpecOption); % Get the file
    [row, ~] = find(DatabaseTimeAltitude(:,1) == TimeInterval); % Find the right row in the file
    if isempty(row)
        row = size(DatabaseAltitudeInc,1); % If the row is empty then take the last row
    end
    if month(startDatum) == month(endDatum) && year(startDatum) == year(endDatum) % Calculate the elapsed time of the month 
        % Start of the month of the current time
        startOfMonth = dateshift(endDatum, 'start', 'month');
        
        % Calculate the elapsed days including hours [hr]
        timeDifference = endDatum - startOfMonth;
        
        daysInCell = days(timeDifference); % [d]
    else
        daysInCell = eomday(year(startDatum), month(startDatum)); % Entire month has elapsed [d]
    end
    CompleteDaysMonth = eomday(year(startDatum), month(startDatum)); % Total number of days in the month

    DateTimeAltitudeInter=DatabaseTimeAltitude(row, 2:end)+((DatabaseTimeAltitude(row+1, 2:end)-DatabaseTimeAltitude(row, 2:end))*daysInCell)/CompleteDaysMonth; % Interpolate month with data of INC = 5°
    
    ScaleMatrix = 1+((DateTimeAltitudeInter(1,:)-DatabaseTimeAltitude(2,2:end))./DatabaseTimeAltitude(2,2:end)); % Determine the percentage by which the data needs to be increased
    
    AktData(:,19) = AktData(:,19) .* ScaleMatrix(1,:)'; % Scaling the value
end

function SurfFlux = ScaleSurfaceFlux(AktData,Surface,SpecOption,i,SMA,INC) % Scaling flux data for the oriented plates
    if INC > 80 && INC < 100
        if strcmp(SpecOption,'Cubesat')
            % Assigment of the flux data for the CubeSats
            SurfFlux.Cubesat.Front{i} = Surface.Cubesat.Front .* AktData;
            SurfFlux.Cubesat.Right{i} = Surface.Cubesat.Right .* AktData;
            SurfFlux.Cubesat.Left{i} = Surface.Cubesat.Left .* AktData;
            SurfFlux.Cubesat.Up{i} = Surface.Cubesat.Up .* AktData;
            SurfFlux.Cubesat.Down{i} = Surface.Cubesat.Down .* AktData; 
            SurfFlux.Cubesat.Rear{i} = Surface.Cubesat.Rear .* AktData;
    
            if SMA > 8125 % Scaling the oriented surfaces with linear equations
                
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
                
                p = [
                    -8.77265000000027e-07	0.00848280254166689
                    3.61225000000009e-07	-0.00259042787500007
                    -3.65885000000012e-07	0.00319357304166677
                    4.22400000000005e-07	-0.00334722866666671
                    -1.25715000000005e-07	0.00115094479166670
                    2.97125000000021e-07	-0.00221041004166684
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
            
                % Coefficient matrix p for SurfFlux.Cubesat.Right
    
                p = [
                    0.000169959925000004	-1.10914061770837
                    -0.000150917410000004	1.55866763275004
                    0.000116537160000002	-0.664248881666687
                    0.000278852810000008	-2.01990421375007
                    0.000200422935000003	-1.35935969079169
                    1.91433500000031e-05	0.162418116083308
                   ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
    
                p = [
                    0.000159993680000003	-1.02504021466669
                    0.000125719895000004	-0.749810904791698
                    9.96559900000020e-05	-0.522368378583350
                    0.000137752280000004	-0.841268778000033
                    0.000430400350000007	-3.25234761558339
                    2.23673500000035e-05	0.138339045416638
                 ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
      
                p = [
                    8.90150000000055e-07	0.0523541137499995
                    6.08239350000015e-05	-0.461070577125012
                    5.16492950000007e-05	-0.390753778125005
                    7.81008300000012e-05	-0.609980621916677
                    4.44674300000004e-05	-0.333256527250003
                    2.49408150000006e-05	-0.1717319214583388
                 ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
       
                p = [
                    2.72373999999995e-06	0.0327505575000004
                    7.46818550000018e-05	-0.577820326458348
                    3.96787750000006e-05	-0.294726415791671
                    8.17323200000014e-05	-0.640466459000011
                    7.40152999999964e-06	-0.0258696310833303
                    3.54348350000009e-05	-0.259626085291674
                 ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
            elseif SMA <= 8125 && SMA > 7025
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
                
                p = [
                    -1.89168727272727e-07	0.00205865874545455
                    -3.78229363636365e-07	0.00328526215681819
                    -9.92429181818183e-08	0.000941600496136365
                    -1.09919236363637e-07	0.000984541397272729
                    -9.19382454545457e-08	0.000809630145681820
                    -1.15596554545455e-07	0.00100480512795455
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Right
                
                p = [
                    -1.36297245454546e-05	0.363121566159091
                    1.32520220000002e-05	0.147116165713635
                    1.52702454545406e-07	0.250302077179546
                    -2.26759290909090e-06	0.267084413195455
                    3.06947036363631e-06	0.236799116086364
                    3.77511890000001e-05	-0.0348494393113644
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
                
                p = [
                    -1.25054615454546e-05	0.353924211752273
                    8.10451281818200e-06	0.191487633493180
                    -8.31779600000001e-06	0.314080306790909
                    1.46160081818188e-06	0.247667539256818
                    1.93140372727269e-06	0.243103695584091
                    2.72067074545457e-05	0.0625805773954530
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
                
                p = [
                    1.53099894545455e-05	-0.0988381733909096
                    1.30664369090910e-05	-0.0854394569500003
                    9.47292200000004e-06	-0.0606238588772730
                    1.21923160000001e-05	-0.0810067422454549
                    9.85266309090911e-06	-0.0657594549136365
                    1.07369470000001e-05	-0.0716653340704549
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
                
                p = [
                    1.44610175454546e-05	-0.0933717534522732
                    1.30143580909091e-05	-0.0851264862659095
                    8.46050654545459e-06	-0.0529965406272730
                    1.16312078181819e-05	-0.0769127107681822
                    1.04432156363637e-05	-0.0700229875363638
                    9.41592572727277e-06	-0.0619645912022731
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
    
             elseif SMA <= 7025 && SMA> 6558
    
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
    
                p = [
                    -5.38965000033973e-11	1.13376396256893e-06	-0.00794909737140370	18.5763225341681
                    -9.62875000059261e-11	1.98374583762024e-06	-0.0136241450648757	31.1924178845907
                    8.56850000047665e-12	-1.61808462509670e-07	0.00100880031600289	-2.07153474507705
                    -1.30123500007938e-10	2.66047058766107e-06	-0.0181297069001518	41.1772144576818
                    -1.80925000011177e-10	3.73985702522678e-06	-0.0257675971359087	59.1774834685039
                    8.03926666716093e-11	-1.65491585010029e-06	0.0113539731952615	-25.9613342431225
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Right
    
                p = [
                    2.79531743350373e-08	-0.000571162816709573	3.88897220420985	-8823.46133224484
                    -2.33795266681509e-09	4.94647385530116e-05	-0.348989114164951	821.340864913194
                    -5.78568600035752e-09	0.000119600499907254	-0.824399311357809	1895.11500174008
                    4.34245833361920e-10	-9.51811626308006e-06	0.0689833771915274	-165.197725448399
                    -2.04725905012608e-08	0.000422010354288082	-2.89980181194269	6642.42602739008
                    -1.71159821677152e-08	0.000351215257158775	-2.40228747722878	5477.46763292226
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
    
                p = [
                    2.90207346684352e-08	-0.000592854871335884	4.03588177462351	-9155.08149769393
                    -3.31387466687459e-09	6.94133384042189e-05	-0.484882070006866	1129.85170846007
                    -5.04392066697901e-09	0.000104453407506337	-0.721301004564942	1661.21365756061
                    5.14174833367671e-10	-1.14042055381968e-05	0.0835827511466918	-202.391321607631
                    -2.08228190012822e-08	0.000429196944451017	-2.94894948523407	6754.44476710433
                    -1.76757991677494e-08	0.000362647321334467	-2.48008857838283	5653.92171918782
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
    
                p = [
                    1.49779116675745e-09	-3.04569007643419e-05	0.206280822374227	-465.308557600069
                    -3.96425450024327e-09	8.14674167674359e-05	-0.558047897168068	1274.17449830192
                    1.99246533345487e-09	-4.07349345524659e-05	0.277525636658343	-630.072761978182
                    1.70509200010457e-09	-3.50239866521217e-05	0.239757898849348	-546.970202336748
                    -1.09634550006818e-09	2.27928185888834e-05	-0.157959748507168	364.919542017628
                    1.03530000004054e-11	-1.46142575008220e-07	0.000534622644430545	-0.125002808140731
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
    
                p = [
                    1.89359533344900e-09	-3.87635800523469e-05	0.264403413148788	-600.900914143493
                    -2.99599466685037e-09	6.15251885537273e-05	-0.421144010676040	960.903929477286
                    2.08859050012738e-09	-4.26942377150845e-05	0.290836272889666	-660.211926177795
                    1.90768266678367e-09	-3.91875005523739e-05	0.268275818051888	-612.071678843906
                    -1.04459750006505e-09	2.17437912388200e-05	-0.150874570096739	348.976206435567
                    1.16405000004643e-11	-1.66781512509414e-07	0.000634723018501118	-0.256036338354232
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
    
            elseif SMA <= 6558 
                
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
                
                p = [
                    0.003585225;
                    0.002484431;
                    0.003173469;
                    0.003285147;
                    0.001745178;
                    0.002002627
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * p(j);
                end
            
                % Coefficient matrix p for SurfFlux.Cubesat.Right
                
                p = [
                    0.503219248;
                    0.504338948;
                    0.48579932;
                    0.458900227;
                    0.468527919;
                    0.486765003
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
                p = [
                    0.503558116;
                    0.504338948;
                    0.48579932;
                    0.458900227;
                    0.468426396;
                    0.490402101
                    ];            
                
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
    
                p = [
                    0.092917655;
                    0.046299132;
                    0.021369048;
                    0.019920635;
                    0.019340102;
                    0.020367751
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
    
                p = [
                    0.036563877;
                    0.027820316;
                    0.022389456;
                    0.020728458;
                    0.020010152;
                    0.02135785
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * p(j,1);
                end
    
            end
        elseif strcmp(SpecOption,'Cubesat with Dragsail')
            % Assigment of the flux data for the CubeSats and DragSail steps
    
            SurfFlux.Cubesat.Front{i} = Surface.Cubesat.Front .* AktData;
            SurfFlux.Cubesat.Right{i} = Surface.Cubesat.Right .* AktData; 
            SurfFlux.Cubesat.Left{i} = Surface.Cubesat.Left .* AktData;
            SurfFlux.Cubesat.Up{i} = Surface.Cubesat.Up .* AktData; 
            SurfFlux.Cubesat.Down{i} = Surface.Cubesat.Down .* AktData; 
            SurfFlux.Cubesat.Rear{i} = Surface.Cubesat.Rear .* AktData; 
            if SMA > 8125 % Linear area
    
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
    
                p = [
                    -8.77265000000027e-07	0.00848280254166689
                    3.61225000000009e-07	-0.00259042787500007
                    -3.65885000000012e-07	0.00319357304166677
                    4.22400000000005e-07	-0.00334722866666671
                    -1.25715000000005e-07	0.00115094479166670
                    2.97125000000021e-07	-0.00221041004166684
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
    
                % Coefficient matrix p for SurfFlux.Cubesat.Right
    
                p = [
                    0.000169959925000004	-1.10914061770837
                    -0.000150917410000004	1.55866763275004
                    0.000116537160000002	-0.664248881666687
                    0.000278852810000008	-2.01990421375007
                    0.000200422935000003	-1.35935969079169
                    1.91433500000031e-05	0.162418116083308
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
                
                p = [
                    0.000159993680000003	-1.02504021466669
                    0.000125719895000004	-0.749810904791698
                    9.96559900000020e-05	-0.522368378583350
                    0.000137752280000004	-0.841268778000033
                    0.000430400350000007	-3.25234761558339
                    2.23673500000035e-05	0.138339045416638
                    ];
                 
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
                
                p = [
                    8.90150000000055e-07	0.0523541137499995
                    6.08239350000015e-05	-0.461070577125012
                    5.16492950000007e-05	-0.390753778125005
                    7.81008300000012e-05	-0.609980621916677
                    4.44674300000004e-05	-0.333256527250003
                    2.49408150000006e-05	-0.1717319214583388
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
    
                p = [
                    2.72373999999995e-06	0.0327505575000004
                    7.46818550000018e-05	-0.577820326458348
                    3.96787750000006e-05	-0.294726415791671
                    8.17323200000014e-05	-0.640466459000011
                    7.40152999999964e-06	-0.0258696310833303
                    3.54348350000009e-05	-0.259626085291674
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
    
            elseif SMA <= 8125 && SMA > 7025 % Linear area
                
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
                
                p = [
                    -1.89168727272727e-07	0.00205865874545455
                    -3.78229363636365e-07	0.00328526215681819
                    -9.92429181818183e-08	0.000941600496136365
                    -1.09919236363637e-07	0.000984541397272729
                    -9.19382454545457e-08	0.000809630145681820
                    -1.15596554545455e-07	0.00100480512795455
                    ];
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Right
                
                p = [
                    -1.36297245454546e-05	0.363121566159091
                    1.32520220000002e-05	0.147116165713635
                    1.52702454545406e-07	0.250302077179546
                    -2.26759290909090e-06	0.267084413195455
                    3.06947036363631e-06	0.236799116086364
                    3.77511890000001e-05	-0.0348494393113644
                    ];
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
                
                p = [
                    -1.25054615454546e-05	0.353924211752273
                    8.10451281818200e-06	0.191487633493180
                    -8.31779600000001e-06	0.314080306790909
                    1.46160081818188e-06	0.247667539256818
                    1.93140372727269e-06	0.243103695584091
                    2.72067074545457e-05	0.0625805773954530
                    ];
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
                
                p = [
                    1.53099894545455e-05	-0.0988381733909096
                    1.30664369090910e-05	-0.0854394569500003
                    9.47292200000004e-06	-0.0606238588772730
                    1.21923160000001e-05	-0.0810067422454549
                    9.85266309090911e-06	-0.0657594549136365
                    1.07369470000001e-05	-0.0716653340704549
                    ];
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
                
                p = [
                    1.44610175454546e-05	-0.0933717534522732
                    1.30143580909091e-05	-0.0851264862659095
                    8.46050654545459e-06	-0.0529965406272730
                    1.16312078181819e-05	-0.0769127107681822
                    1.04432156363637e-05	-0.0700229875363638
                    9.41592572727277e-06	-0.0619645912022731
                    ];
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
    
             elseif SMA <= 7025 && SMA> 6558 %Polynomial area
    
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
    
                p = [
                    -5.38965000033973e-11	1.13376396256893e-06	-0.00794909737140370	18.5763225341681
                    -9.62875000059261e-11	1.98374583762024e-06	-0.0136241450648757	31.1924178845907
                    8.56850000047665e-12	-1.61808462509670e-07	0.00100880031600289	-2.07153474507705
                    -1.30123500007938e-10	2.66047058766107e-06	-0.0181297069001518	41.1772144576818
                    -1.80925000011177e-10	3.73985702522678e-06	-0.0257675971359087	59.1774834685039
                    8.03926666716093e-11	-1.65491585010029e-06	0.0113539731952615	-25.9613342431225
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Right
    
                p = [
                    2.79531743350373e-08	-0.000571162816709573	3.88897220420985	-8823.46133224484
                    -2.33795266681509e-09	4.94647385530116e-05	-0.348989114164951	821.340864913194
                    -5.78568600035752e-09	0.000119600499907254	-0.824399311357809	1895.11500174008
                    4.34245833361920e-10	-9.51811626308006e-06	0.0689833771915274	-165.197725448399
                    -2.04725905012608e-08	0.000422010354288082	-2.89980181194269	6642.42602739008
                    -1.71159821677152e-08	0.000351215257158775	-2.40228747722878	5477.46763292226
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
    
                p = [
                    2.90207346684352e-08	-0.000592854871335884	4.03588177462351	-9155.08149769393
                    -3.31387466687459e-09	6.94133384042189e-05	-0.484882070006866	1129.85170846007
                    -5.04392066697901e-09	0.000104453407506337	-0.721301004564942	1661.21365756061
                    5.14174833367671e-10	-1.14042055381968e-05	0.0835827511466918	-202.391321607631
                    -2.08228190012822e-08	0.000429196944451017	-2.94894948523407	6754.44476710433
                    -1.76757991677494e-08	0.000362647321334467	-2.48008857838283	5653.92171918782
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
    
                p = [
                    1.49779116675745e-09	-3.04569007643419e-05	0.206280822374227	-465.308557600069
                    -3.96425450024327e-09	8.14674167674359e-05	-0.558047897168068	1274.17449830192
                    1.99246533345487e-09	-4.07349345524659e-05	0.277525636658343	-630.072761978182
                    1.70509200010457e-09	-3.50239866521217e-05	0.239757898849348	-546.970202336748
                    -1.09634550006818e-09	2.27928185888834e-05	-0.157959748507168	364.919542017628
                    1.03530000004054e-11	-1.46142575008220e-07	0.000534622644430545	-0.125002808140731
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
    
                p = [
                    1.89359533344900e-09	-3.87635800523469e-05	0.264403413148788	-600.900914143493
                    -2.99599466685037e-09	6.15251885537273e-05	-0.421144010676040	960.903929477286
                    2.08859050012738e-09	-4.26942377150845e-05	0.290836272889666	-660.211926177795
                    1.90768266678367e-09	-3.91875005523739e-05	0.268275818051888	-612.071678843906
                    -1.04459750006505e-09	2.17437912388200e-05	-0.150874570096739	348.976206435567
                    1.16405000004643e-11	-1.66781512509414e-07	0.000634723018501118	-0.256036338354232
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
    
            elseif SMA <= 6558 % Constant area
    
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
    
                p = [
                    0.003585225;
                    0.002484431;
                    0.003173469;
                    0.003285147;
                    0.001745178;
                    0.002002627
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Right
    
                p = [
                    0.503219248;
                    0.504338948;
                    0.48579932;
                    0.458900227;
                    0.468527919;
                    0.486765003
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
    
                p = [
                    0.503558116;
                    0.504338948;
                    0.48579932;
                    0.458900227;
                    0.468426396;
                    0.490402101
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
    
                p = [
                    0.092917655;
                    0.046299132;
                    0.021369048;
                    0.019920635;
                    0.019340102;
                    0.020367751
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
    
                p = [
                    0.036563877;
                    0.027820316;
                    0.022389456;
                    0.020728458;
                    0.020010152;
                    0.02135785
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * p(j,1);
                end
    
            end
            if i > 1
                SurfFlux.Dragsail.Front{i-1} = Surface.Dragsail.Front{i-1} .* AktData;
                SurfFlux.Dragsail.Rear{i-1} = Surface.Dragsail.Rear{i-1} .* AktData; 
                if SMA > 8125 % Scaling the oriented surfaces with linear equations
    
                    % Coefficient matrix p for SurfFlux.Dragsail.Rear
    
                    p = [
                        -8.77265000000027e-07	0.00848280254166689
                        3.61225000000009e-07	-0.00259042787500007
                        -3.65885000000012e-07	0.00319357304166677
                        4.22400000000005e-07	-0.00334722866666671
                        -1.25715000000005e-07	0.00115094479166670
                        2.97125000000021e-07	-0.00221041004166684
                        ];
    
                    for j = 1:6
                        SurfFlux.Dragsail.Rear{i-1}(j,19) = SurfFlux.Dragsail.Rear{i-1}(j,19) * (p(j,1)*SMA + p(j,2));
                    end
    
               elseif SMA <= 8125 && SMA > 7025
    
                    % Coefficient matrix p for SurfFlux.Dragsail.Rear
    
                    p = [
                        -1.89168727272727e-07	0.00205865874545455
                        -3.78229363636365e-07	0.00328526215681819
                        -9.92429181818183e-08	0.000941600496136365
                        -1.09919236363637e-07	0.000984541397272729
                        -9.19382454545457e-08	0.000809630145681820
                        -1.15596554545455e-07	0.00100480512795455
                        ];
    
                    for j = 1:6
                        SurfFlux.Dragsail.Rear{i-1}(j,19) = SurfFlux.Dragsail.Rear{i-1}(j,19) * (p(j,1)*SMA + p(j,2));
                    end
    
                elseif SMA <= 7025 && SMA> 6558
    
                    % Coefficient matrix p for SurfFlux.Dragsail.Rear
    
                    p = [
                        -5.38965000033973e-11	1.13376396256893e-06	-0.00794909737140370	18.5763225341681
                        -9.62875000059261e-11	1.98374583762024e-06	-0.0136241450648757	31.1924178845907
                        8.56850000047665e-12	-1.61808462509670e-07	0.00100880031600289	-2.07153474507705
                        -1.30123500007938e-10	2.66047058766107e-06	-0.0181297069001518	41.1772144576818
                        -1.80925000011177e-10	3.73985702522678e-06	-0.0257675971359087	59.1774834685039
                        8.03926666716093e-11	-1.65491585010029e-06	0.0113539731952615	-25.9613342431225
                        ];
    
                    for j = 1:6
                        SurfFlux.Dragsail.Rear{i-1}(j,19) = SurfFlux.Dragsail.Rear{i-1}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                    end
    
                elseif SMA <= 6558   
    
                    % Coefficient matrix p for SurfFlux.Dragsail.Rear
    
                    p = [
                        0.003585225;
                        0.002484431;
                        0.003173469;
                        0.003285147;
                        0.001745178;
                        0.002002627
                        ];
    
                    for j = 1:6
                        SurfFlux.Dragsail.Rear{i-1}(j,19) = SurfFlux.Dragsail.Rear{i-1}(j,19) * p(j,1);
                    end
    
                end
            end
        else
            SurfFlux.Sphere{i} = Surface.Sphere{1} .* (AktData./4);
        end
    else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(SpecOption,'Cubesat')
            % Zuweisung des Fluxes fuer die Cubesats
            SurfFlux.Cubesat.Front{i} = Surface.Cubesat.Front .* AktData;
            SurfFlux.Cubesat.Right{i} = Surface.Cubesat.Right .* AktData;
            SurfFlux.Cubesat.Left{i} = Surface.Cubesat.Left .* AktData;
            SurfFlux.Cubesat.Up{i} = Surface.Cubesat.Up .* AktData;
            SurfFlux.Cubesat.Down{i} = Surface.Cubesat.Down .* AktData; 
            SurfFlux.Cubesat.Rear{i} = Surface.Cubesat.Rear .* AktData;
    
            if SMA > 8125 % Scaling the oriented surfaces with linear equations
                
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
                
                p = [
                    -2.45e-5,  2.644e-1;
                    3.892e-5, -2.706e-1;
                    8.293e-5, -6.261e-1;
                    3.396e-4, -2.752;
                    1.515e-4, -1.21;
                    5.44e-5,  -3.965e-1
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
            
                % Coefficient matrix p for SurfFlux.Cubesat.Right
    
                p = [
                    6.404e-4,  -4.594;
                    8.452e-4,  -6.328;
                    4.975e-4,  -3.422;
                   -1.796e-4,   2.161;
                   -2.792e-4,   3.054;
                   -3.979e-4,   4.071
                   ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
    
                p = [
                -1.158e-3,  1.028e1;
                -3.583e-4,  3.662;
                 4.871e-4, -3.25;
                 1.321e-3, -1.018e1;
                 1.027e-3, -7.891;
                 2.43e-4,   -1.386
                 ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
      
                p = [
                -7.077e-5,  7.819e-1;
                 3.104e-5, -1.294e-1;
                 2.108e-4, -1.597;
                 4.134e-4, -3.305;
                 2.843e-4, -2.236;
                 1.962e-4, -1.488
                 ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
       
                p = [
                 1.636e-6,  1.573e-1;
                -3.151e-5,  3.953e-1;
                 1.279e-4, -9.249e-1;
                 2.974e-4, -2.323;
                 1.532e-4, -1.168;
                 1.233e-4, -9.026e-1
                 ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
            elseif SMA <= 8125 && SMA > 7025
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
                
                p = [
                     6.129e-6, -2.673e-2;
                     1.059e-5, -6.569e-2;
                     5.020e-6, -2.241e-2;
                     1.219e-5, -7.990e-2;
                     5.298e-6, -3.325e-2;
                     1.405e-5, -9.118e-2;
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Right
                
                p = [
                    -5.712e-5,  9.935e-1;
                     3.868e-5,  2.631e-1;
                     8.166e-6,  4.670e-1;
                     2.989e-5,  2.874e-1;
                    -1.204e-5,  5.927e-1;
                     4.938e-5,  2.035e-1;
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
                
                p = [
                     5.380e-5,  1.747e-1;
                     4.885e-5,  1.931e-1;
                     3.269e-5,  3.565e-1;
                     1.608e-4, -6.277e-1;
                     6.401e-5,  7.899e-2;
                     1.505e-4, -6.111e-1;
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
                
                p = [
                     4.626e-5, -3.062e-1;
                     3.690e-5, -2.479e-1;
                     2.602e-5, -1.679e-1;
                     3.441e-5, -2.343e-1;
                     2.047e-5, -1.361e-1;
                     3.348e-5, -2.280e-1;
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
                
                p = [
                     3.713e-5, -2.416e-1;
                     3.868e-5, -2.631e-1;
                     3.131e-5, -2.066e-1;
                     3.417e-5, -2.321e-1;
                     1.906e-5, -1.264e-1;
                     4.022e-5, -2.751e-1;
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
    
             elseif SMA <= 7025 && SMA> 6558
    
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
    
                p = [-7.32559945007295e-08,0.00150496076678445,-10.3063884909542,23528.3362872214;
                    -1.63677412667677e-07,0.00336525300232765,-23.0621767880850,52679.0214518787;
                    -4.11602179919497e-07,0.00845259585668511,-57.8551058196769,131987.564527940;
                    -1.06721638000839e-07,0.00218983902643132,-14.9767727815125,34140.6631292157;
                    7.27563929174125e-08,-0.00147356904480353,9.94247872480426,-22347.7183119919;
                    -5.31851231420567e-07,0.0109032793895086,-74.4982997462064,169651.609514051
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Right
    
                p = [-1.51507613584882e-07,0.00311377954331269,-21.3328154241673,48721.8426198415;
                    -2.18183804997732e-08,0.000459115684336458,-3.22197136254680,7541.31547829916;
                    -4.50107945836229e-08,0.000938430654197841,-6.52208258412019,15110.7271868036;
                    1.79901140835078e-08,-0.000356280004618719,2.34658827593987,-5138.30281564287;
                    -5.87312656669225e-08,0.00120316909389447,-8.21804187371506,18716.0545534252;
                    -1.49034404667662e-07,0.00306248052034162,-20.9751135147865,47883.3885669635
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
    
                p = [-2.32705589171277e-08,0.000473662639240615,-3.21542051452970,7280.74098494595;
                    2.20881508334584e-08,-0.000435057586204321,2.84824751874017,-6195.36490158282;
                    -1.14585435501143e-07,0.00234946212688213,-16.0575475572739,36582.3407211498;
                    -1.08796366084316e-07,0.00221857008059226,-15.0781479669276,34154.4625428346;
                    2.49428323919222e-07,-0.00507386206366526,34.3933298013259,-77686.6656260824;
                    -3.09793913085669e-07,0.00634727143548043,-43.3440218113853,98651.2724564319
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
    
                p = [-2.09704937502633e-08,0.000430601363609812,-2.94785668881340,6728.38147640328;
                    -2.32277525834521e-08,0.000479350620522945,-3.29775106855762,7563.19528176778;
                    -3.83623071669234e-08,0.000789510950374854,-5.41604261473129,12384.5025630829;
                    -1.54754827501411e-08,0.000317337705168042,-2.16906727982488,4941.99799163881;
                    3.54085286670134e-08,-0.000719924369328466,4.87715003478505,-11008.7947235688;
                    -5.71075250003949e-08,0.00117089735627587,-8.00153106946428,18224.5622153244
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
    
                p = [-1.04739225001151e-08,0.000215095483393409,-1.47277891260015,3362.32609874467;
                    -2.12668423334036e-08,0.000438343814183566,-3.01169112625494,6897.52668051556;
                    -5.29734823337121e-08,0.00108782717618983,-7.44586023970730,16987.3350034181;
                    -7.88211966673720e-09,0.000161943171412146,-1.10917849196382,2532.59763334825;
                    1.77493031668469e-08,-0.000360055120573301,2.43333070402330,-5478.57642173795;
                    -6.66978004171308e-08,0.00136690870734245,-9.33668480389638,21255.5286409980
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
    
            elseif SMA <= 6558 
                
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
                
                p = [
                    1.8320;
                    4.2608;
                    3.2551;
                    3.5442;
                    0.9473;
                    8.3434
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * p(j);
                end
            
                % Coefficient matrix p for SurfFlux.Cubesat.Right
                
                p = [
                    4.4536;
                    3.8499;
                    1.4325;
                    3.5966;
                    1.9683;
                    3.9515
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
                p = [
                    1.6189;
                    3.0972;
                    1.3140;
                    2.7686;
                    0.8524;
                    5.6733
                    ];            
                
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
    
                p = [
                    0.7051;
                    1.0529;
                    0.4338;
                    0.6268;
                    0.2059;
                    0.9844
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
    
                p = [
                    0.3915;
                    0.7899;
                    0.4720;
                    0.4530;
                    0.2074;
                    1.0727
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * p(j,1);
                end
    
            end
        elseif strcmp(SpecOption,'Cubesat with Dragsail')
            % Assigment of the flux data for the CubeSats and DragSail steps
    
            SurfFlux.Cubesat.Front{i} = Surface.Cubesat.Front .* AktData;
            SurfFlux.Cubesat.Right{i} = Surface.Cubesat.Right .* AktData; 
            SurfFlux.Cubesat.Left{i} = Surface.Cubesat.Left .* AktData;
            SurfFlux.Cubesat.Up{i} = Surface.Cubesat.Up .* AktData; 
            SurfFlux.Cubesat.Down{i} = Surface.Cubesat.Down .* AktData; 
            SurfFlux.Cubesat.Rear{i} = Surface.Cubesat.Rear .* AktData; 
            if SMA > 8125 % Linear area
    
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
    
                p = [
                    -2.45e-5,  2.644e-1;
                     3.892e-5, -2.706e-1;
                     8.293e-5, -6.261e-1;
                     3.396e-4, -2.752;
                     1.515e-4, -1.21;
                     5.44e-5,  -3.965e-1
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
    
                % Coefficient matrix p for SurfFlux.Cubesat.Right
    
                p = [
                     6.404e-4,  -4.594;
                     8.452e-4,  -6.328;
                     4.975e-4,  -3.422;
                    -1.796e-4,   2.161;
                    -2.792e-4,   3.054;
                    -3.979e-4,   4.071
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
                
                p = [
                    -1.158e-3,  10.28;
                    -3.583e-4,   3.662;
                     4.871e-4,  -3.25;
                     1.321e-3, -10.18;
                     1.027e-3,  -7.891;
                     2.43e-4,   -1.386
                    ];
                 
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
                
                p = [
                    -7.077e-5,  7.819e-1;
                     3.104e-5, -1.294e-1;
                     2.108e-4, -1.597;
                     4.134e-4, -3.305;
                     2.843e-4, -2.236;
                     1.962e-4, -1.488
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
    
                p = [
                     1.636e-6,  1.573e-1;
                    -3.151e-5,  3.953e-1;
                     1.279e-4, -9.249e-1;
                     2.974e-4, -2.323;
                     1.532e-4, -1.168;
                     1.233e-4, -9.026e-1
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
    
            elseif SMA <= 8125 && SMA > 7025 % Linear area
                
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
                
                p = [
                     6.129e-6, -2.673e-2;
                     1.059e-5, -6.569e-2;
                     5.020e-6, -2.241e-2;
                     1.219e-5, -7.990e-2;
                     5.298e-6, -3.325e-2;
                     1.405e-5, -9.118e-2
                    ];
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Right
                
                p = [
                    -5.712e-5,  9.935e-1;
                     3.868e-5,  2.631e-1;
                     8.166e-6,  4.670e-1;
                     2.989e-5,  2.874e-1;
                    -1.204e-5,  5.927e-1;
                     4.938e-5,  2.035e-1
                    ];
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
                
                p = [
                     5.380e-5,  1.747e-1;
                     4.885e-5,  1.931e-1;
                     3.269e-5,  3.565e-1;
                     1.608e-4, -6.277e-1;
                     6.401e-5,  7.899e-2;
                     1.505e-4, -6.111e-1
                    ];
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
                
                p = [
                     4.626e-5, -3.062e-1;
                     3.690e-5, -2.479e-1;
                     2.602e-5, -1.679e-1;
                     3.441e-5, -2.343e-1;
                     2.047e-5, -1.361e-1;
                     3.348e-5, -2.280e-1
                    ];
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
                
                p = [
                     3.713e-5, -2.416e-1;
                     3.868e-5, -2.631e-1;
                     3.131e-5, -2.066e-1;
                     3.417e-5, -2.321e-1;
                     1.906e-5, -1.264e-1;
                     4.022e-5, -2.751e-1
                    ];
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1)*SMA + p(j,2));
                end
    
             elseif SMA <= 7025 && SMA> 6558 %Polynomial area
    
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
    
                p = [-7.32559945007295e-08,0.00150496076678445,-10.3063884909542,23528.3362872214;
                    -1.63677412667677e-07,0.00336525300232765,-23.0621767880850,52679.0214518787;
                    -4.11602179919497e-07,0.00845259585668511,-57.8551058196769,131987.564527940;
                    -1.06721638000839e-07,0.00218983902643132,-14.9767727815125,34140.6631292157;
                    7.27563929174125e-08,-0.00147356904480353,9.94247872480426,-22347.7183119919;
                    -5.31851231420567e-07,0.0109032793895086,-74.4982997462064,169651.609514051
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Right
    
                p = [-1.51507613584882e-07,0.00311377954331269,-21.3328154241673,48721.8426198415;
                    -2.18183804997732e-08,0.000459115684336458,-3.22197136254680,7541.31547829916;
                    -4.50107945836229e-08,0.000938430654197841,-6.52208258412019,15110.7271868036;
                    1.79901140835078e-08,-0.000356280004618719,2.34658827593987,-5138.30281564287;
                    -5.87312656669225e-08,0.00120316909389447,-8.21804187371506,18716.0545534252;
                    -1.49034404667662e-07,0.00306248052034162,-20.9751135147865,47883.3885669635
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
    
                p = [-2.32705589171277e-08,0.000473662639240615,-3.21542051452970,7280.74098494595;
                    2.20881508334584e-08,-0.000435057586204321,2.84824751874017,-6195.36490158282;
                    -1.14585435501143e-07,0.00234946212688213,-16.0575475572739,36582.3407211498;
                    -1.08796366084316e-07,0.00221857008059226,-15.0781479669276,34154.4625428346;
                    2.49428323919222e-07,-0.00507386206366526,34.3933298013259,-77686.6656260824;
                    -3.09793913085669e-07,0.00634727143548043,-43.3440218113853,98651.2724564319
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
    
                p = [-2.09704937502633e-08,0.000430601363609812,-2.94785668881340,6728.38147640328;
                    -2.32277525834521e-08,0.000479350620522945,-3.29775106855762,7563.19528176778;
                    -3.83623071669234e-08,0.000789510950374854,-5.41604261473129,12384.5025630829;
                    -1.54754827501411e-08,0.000317337705168042,-2.16906727982488,4941.99799163881;
                    3.54085286670134e-08,-0.000719924369328466,4.87715003478505,-11008.7947235688;
                    -5.71075250003949e-08,0.00117089735627587,-8.00153106946428,18224.5622153244
                    ];
                
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
    
                p = [-1.04739225001151e-08,0.000215095483393409,-1.47277891260015,3362.32609874467;
                    -2.12668423334036e-08,0.000438343814183566,-3.01169112625494,6897.52668051556;
                    -5.29734823337121e-08,0.00108782717618983,-7.44586023970730,16987.3350034181;
                    -7.88211966673720e-09,0.000161943171412146,-1.10917849196382,2532.59763334825;
                    1.77493031668469e-08,-0.000360055120573301,2.43333070402330,-5478.57642173795;
                    -6.66978004171308e-08,0.00136690870734245,-9.33668480389638,21255.5286409980];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                end
    
            elseif SMA <= 6558 % Constant area
    
                % Coefficient matrix p for SurfFlux.Cubesat.Rear
    
                p = [
                    1.8320;
                    4.2608;
                    3.2551;
                    3.5442;
                    0.9473;
                    8.3434
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Rear{i}(j,19) = SurfFlux.Cubesat.Rear{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Right
    
                p = [
                    4.4536;
                    3.8499;
                    1.4325;
                    3.5966;
                    1.9683;
                    3.9515
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Right{i}(j,19) = SurfFlux.Cubesat.Right{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Left
    
                p = [
                    1.6189;
                    3.0972;
                    1.3140;
                    2.7686;
                    0.8524;
                    5.6733
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Left{i}(j,19) = SurfFlux.Cubesat.Left{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Up
    
                p = [
                    0.7051;
                    1.0529;
                    0.4338;
                    0.6268;
                    0.2059;
                    0.9844
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Up{i}(j,19) = SurfFlux.Cubesat.Up{i}(j,19) * p(j,1);
                end
                
                % Coefficient matrix p for SurfFlux.Cubesat.Down
    
                p = [
                    0.3915;
                    0.7899;
                    0.4720;
                    0.4530;
                    0.2074;
                    1.0727
                    ];
    
                for j = 1:6
                    SurfFlux.Cubesat.Down{i}(j,19) = SurfFlux.Cubesat.Down{i}(j,19) * p(j,1);
                end
    
            end
            if i > 1
                SurfFlux.Dragsail.Front{i-1} = Surface.Dragsail.Front{i-1} .* AktData;
                SurfFlux.Dragsail.Rear{i-1} = Surface.Dragsail.Rear{i-1} .* AktData; 
                if SMA > 8125 % Scaling the oriented surfaces with linear equations
    
                    % Coefficient matrix p for SurfFlux.Dragsail.Rear
    
                    p = [
                        -2.45e-5,  2.644e-1;
                         3.892e-5, -2.706e-1;
                         8.293e-5, -6.261e-1;
                         3.396e-4, -2.752;
                         1.515e-4, -1.21;
                         5.44e-5,  -3.965e-1
                        ];
    
                    for j = 1:6
                        SurfFlux.Dragsail.Rear{i-1}(j,19) = SurfFlux.Dragsail.Rear{i-1}(j,19) * (p(j,1)*SMA + p(j,2));
                    end
    
               elseif SMA <= 8125 && SMA > 7025
    
                    % Coefficient matrix p for SurfFlux.Dragsail.Rear
    
                    p = [
                        6.129e-6,  -2.673e-2;
                        1.059e-5,  -6.569e-2;
                        5.020e-6,  -2.241e-2;
                        1.219e-5,  -7.990e-2;
                        5.298e-6,  -3.325e-2;
                        1.405e-5,  -9.118e-2
                        ];
    
                    for j = 1:6
                        SurfFlux.Dragsail.Rear{i-1}(j,19) = SurfFlux.Dragsail.Rear{i-1}(j,19) * (p(j,1)*SMA + p(j,2));
                    end
    
                elseif SMA <= 7025 && SMA> 6558
    
                    % Coefficient matrix p for SurfFlux.Dragsail.Rear
    
                    p = [-7.32559945007295e-08,0.00150496076678445,-10.3063884909542,23528.3362872214;
                        -1.63677412667677e-07,0.00336525300232765,-23.0621767880850,52679.0214518787;
                        -4.11602179919497e-07,0.00845259585668511,-57.8551058196769,131987.564527940;
                        -1.06721638000839e-07,0.00218983902643132,-14.9767727815125,34140.6631292157;
                        7.27563929174125e-08,-0.00147356904480353,9.94247872480426,-22347.7183119919;
                        -5.31851231420567e-07,0.0109032793895086,-74.4982997462064,169651.609514051
                        ];
    
                    for j = 1:6
                        SurfFlux.Dragsail.Rear{i-1}(j,19) = SurfFlux.Dragsail.Rear{i-1}(j,19) * (p(j,1) * SMA^3 + p(j,2) * SMA^2 + p(j,3) * SMA + p(j,4));
                    end
    
                elseif SMA <= 6558   
    
                    % Coefficient matrix p for SurfFlux.Dragsail.Rear
    
                    p = [
                        1.8320;
                        4.2608;
                        3.2551;
                        3.5442;
                        0.9473;
                        8.3434
                        ];
    
                    for j = 1:6
                        SurfFlux.Dragsail.Rear{i-1}(j,19) = SurfFlux.Dragsail.Rear{i-1}(j,19) * p(j,1);
                    end
    
                end
            end
        else
            SurfFlux.Sphere{i} = Surface.Sphere{1} .* (AktData./4);
        end
    end
end

function SurfFlux = ScaleTimeInterval(SurfFlux,t,SpecOption,i) % Scaling flux data for time spent in cell
    if strcmp(SpecOption,'Cubesat')
        % Assigment of the flux data for the CubeSats
        SurfFlux.Cubesat.Front{i} = SurfFlux.Cubesat.Front{i}* t;
        SurfFlux.Cubesat.Right{i} = SurfFlux.Cubesat.Right{i} * t;
        SurfFlux.Cubesat.Left{i} = SurfFlux.Cubesat.Left{i} * t;
        SurfFlux.Cubesat.Up{i} = SurfFlux.Cubesat.Up{i} * t;
        SurfFlux.Cubesat.Down{i} = SurfFlux.Cubesat.Down{i} * t;
        SurfFlux.Cubesat.Rear{i} = SurfFlux.Cubesat.Rear{i} * t;
    elseif strcmp(SpecOption,'Cubesat with Dragsail')
        % Assigment of the flux data for the CubeSats and DragSail steps

        SurfFlux.Cubesat.Front{i} = SurfFlux.Cubesat.Front{i} * t;
        SurfFlux.Cubesat.Right{i} = SurfFlux.Cubesat.Right{i} * t;
        SurfFlux.Cubesat.Left{i} = SurfFlux.Cubesat.Left{i} * t;
        SurfFlux.Cubesat.Up{i} = SurfFlux.Cubesat.Up{i} * t;
        SurfFlux.Cubesat.Down{i} = SurfFlux.Cubesat.Down{i} * t;
        SurfFlux.Cubesat.Rear{i} = SurfFlux.Cubesat.Rear{i} * t;
        if i > 1
            SurfFlux.Dragsail.Front{i-1} = SurfFlux.Dragsail.Front{i-1} * t;
            SurfFlux.Dragsail.Rear{i-1} = SurfFlux.Dragsail.Rear{i-1} * t;
        end
    else
        SurfFlux.Sphere{i} = SurfFlux.Sphere{i} .* t;
    end

end

function SurfResult = AdaptResultSurfFlux(SurfResult, SurfFlux,SpecOption,i) % Adds total flux of previous cells to current cell
    if strcmp(SpecOption,'Cubesat')
        % Assigment of the flux data for the CubeSats
        SurfResult.Cubesat.Front{i} = SurfFlux.Cubesat.Front{i}+ SurfResult.Cubesat.Front{i};
        SurfResult.Cubesat.Right{i} = SurfFlux.Cubesat.Right{i} + SurfResult.Cubesat.Right{i};
        SurfResult.Cubesat.Left{i} = SurfFlux.Cubesat.Left{i} + SurfResult.Cubesat.Left{i};
        SurfResult.Cubesat.Up{i} = SurfFlux.Cubesat.Up{i} + SurfResult.Cubesat.Up{i};
        SurfResult.Cubesat.Down{i} = SurfFlux.Cubesat.Down{i} + SurfResult.Cubesat.Down{i};
        SurfResult.Cubesat.Rear{i} = SurfFlux.Cubesat.Rear{i} + SurfResult.Cubesat.Rear{i};
    elseif strcmp(SpecOption,'Cubesat with Dragsail')
        % Assigment of the flux data for the CubeSats and DragSail steps

        SurfResult.Cubesat.Front{i} = SurfFlux.Cubesat.Front{i} + SurfResult.Cubesat.Front{i};
        SurfResult.Cubesat.Right{i} = SurfFlux.Cubesat.Right{i} + SurfResult.Cubesat.Right{i};
        SurfResult.Cubesat.Left{i} = SurfFlux.Cubesat.Left{i} + SurfResult.Cubesat.Left{i};
        SurfResult.Cubesat.Up{i} = SurfFlux.Cubesat.Up{i} + SurfResult.Cubesat.Up{i};
        SurfResult.Cubesat.Down{i} = SurfFlux.Cubesat.Down{i} + SurfResult.Cubesat.Down{i};
        SurfResult.Cubesat.Rear{i} = SurfFlux.Cubesat.Rear{i} + SurfResult.Cubesat.Rear{i};
        if i > 1
            SurfResult.Dragsail.Front{i-1} = SurfFlux.Dragsail.Front{i-1} + SurfResult.Dragsail.Front{i-1};
            SurfResult.Dragsail.Rear{i-1} = SurfFlux.Dragsail.Rear{i-1} + SurfResult.Dragsail.Rear{i-1};
        end
    else
        SurfResult.Sphere{i} = SurfFlux.Sphere{i} + SurfResult.Sphere{i};
    end
end

function SurfResult = UnitCorrection(SurfResult,SpecOption,nDaysEnd,i) % Corrects the unit back to years
    if strcmp(SpecOption,'Cubesat')
        % Assigment of the flux data for the CubeSats
        SurfResult.Cubesat.Front{i}(:, 2:end) = SurfResult.Cubesat.Front{i}(:, 2:end) ./ (nDaysEnd/365);
        SurfResult.Cubesat.Right{i}(:, 2:end) = SurfResult.Cubesat.Right{i}(:, 2:end) ./ (nDaysEnd/365);
        SurfResult.Cubesat.Left{i}(:, 2:end) = SurfResult.Cubesat.Left{i}(:, 2:end) ./ (nDaysEnd/365);
        SurfResult.Cubesat.Up{i}(:, 2:end) = SurfResult.Cubesat.Up{i}(:, 2:end) ./ (nDaysEnd/365);
        SurfResult.Cubesat.Down{i}(:, 2:end) = SurfResult.Cubesat.Down{i}(:, 2:end) ./ (nDaysEnd/365);
        SurfResult.Cubesat.Rear{i}(:, 2:end) = SurfResult.Cubesat.Rear{i}(:, 2:end) ./ (nDaysEnd/365);
    elseif strcmp(SpecOption,'Cubesat with Dragsail')
        % Assigment of the flux data for the CubeSats and DragSail steps

        SurfResult.Cubesat.Front{i}(:, 2:end) = SurfResult.Cubesat.Front{i}(:, 2:end) ./ (nDaysEnd/365);
        SurfResult.Cubesat.Right{i}(:, 2:end) = SurfResult.Cubesat.Right{i}(:, 2:end) ./ (nDaysEnd/365);
        SurfResult.Cubesat.Left{i}(:, 2:end) = SurfResult.Cubesat.Left{i}(:, 2:end) ./ (nDaysEnd/365);
        SurfResult.Cubesat.Up{i}(:, 2:end) = SurfResult.Cubesat.Up{i}(:, 2:end) ./ (nDaysEnd/365);
        SurfResult.Cubesat.Down{i}(:, 2:end) = SurfResult.Cubesat.Down{i}(:, 2:end) ./ (nDaysEnd/365);
        SurfResult.Cubesat.Rear{i}(:, 2:end) = SurfResult.Cubesat.Rear{i}(:, 2:end) ./ (nDaysEnd/365);
        if i > 1
            SurfResult.Dragsail.Front{i-1}(:, 2:end) = SurfResult.Dragsail.Front{i-1}(:, 2:end) ./ (nDaysEnd/365);
            SurfResult.Dragsail.Rear{i-1}(:, 2:end) = SurfResult.Dragsail.Rear{i-1}(:, 2:end) ./ (nDaysEnd/365);
        end
    else
        SurfResult.Sphere{i}(:, 2:end) = SurfResult.Sphere{i}(:, 2:end) ./ (nDaysEnd/365);
    end
end

function SurfResult = CreateSurfResultStruct(SpecOption,i) % Creates a structure for the SurfResult
    if strcmp(SpecOption,'Cubesat')
        % Assigment of the flux data for the CubeSats
        SurfResult.Cubesat.Front{i} = zeros(6,19); 
        SurfResult.Cubesat.Front{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
        SurfResult.Cubesat.Right{i} = zeros(6,19); 
        SurfResult.Cubesat.Right{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];        
        SurfResult.Cubesat.Left{i} = zeros(6,19); 
        SurfResult.Cubesat.Left{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];        
        SurfResult.Cubesat.Up{i} = zeros(6,19); 
        SurfResult.Cubesat.Up{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];        
        SurfResult.Cubesat.Down{i} = zeros(6,19); 
        SurfResult.Cubesat.Down{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];        
        SurfResult.Cubesat.Rear{i} = zeros(6,19); 
        SurfResult.Cubesat.Rear{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
    elseif strcmp(SpecOption,'Cubesat with Dragsail')
        % Assigment of the flux data for the CubeSats and DragSail steps

        SurfResult.Cubesat.Front{i} = zeros(6,19); 
        SurfResult.Cubesat.Front{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
        SurfResult.Cubesat.Right{i} = zeros(6,19); 
        SurfResult.Cubesat.Right{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];        
        SurfResult.Cubesat.Left{i} = zeros(6,19); 
        SurfResult.Cubesat.Left{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];        
        SurfResult.Cubesat.Up{i} = zeros(6,19); 
        SurfResult.Cubesat.Up{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];        
        SurfResult.Cubesat.Down{i} = zeros(6,19); 
        SurfResult.Cubesat.Down{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];        
        SurfResult.Cubesat.Rear{i} = zeros(6,19); 
        SurfResult.Cubesat.Rear{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];
        if i > 1
            SurfResult.Dragsail.Front{i-1} = zeros(6,19); 
            SurfResult.Dragsail.Front{i-1}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];     
            SurfResult.Dragsail.Rear{i-1} = zeros(6,19);
            SurfResult.Dragsail.Rear{i-1}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];     
        end
    else
        SurfResult.Sphere{i} = zeros(6,19);
        SurfResult.Sphere{i}(:, 1) = [0.003162; 0.03162; 0.3162; 3.162;31.62;316.2];     
    end
end
