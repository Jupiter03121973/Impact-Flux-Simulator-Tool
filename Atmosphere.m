function [rho,varargout] = Atmosphere(AtmospericData,r,Phi,Theta,t,R,varargin)              
    
    persistent rhoV
    persistent tV
    persistent hV
    
    if nargin == 6
        
    atmosphere = AtmospericData.atmosphere;

    if contains(atmosphere, "Exponential")
        %% 
        rho0 = 1;
        altitude = r-R;                     % h in [m] 
        S = [-0.160 -0.160 -0.110 -0.085 -0.067 -0.057 -0.050 -0.045 -0.040 -0.037 -0.035]'.*1e-3; % Scale heights [1/km]
        h = [0 100 200 300 400 500 600 700 800 900 1000]'.*1000; % Scale height borders                

        % Determine scale height
        switch true
            case altitude < h(1)
                rho = rho0*exp(S(1)*altitude);
            case h(1) <= altitude && altitude < h(1+1)
                Seff = S(1)+((altitude-h(1))/100000)*(S(1+1)-S(1));
                rho = rho0*exp(Seff*altitude);
            case h(2) <= altitude && altitude < h(2+1)
                Seff = S(2)+((altitude-h(2))/100000)*(S(2+1)-S(2));
                rho = rho0*exp(Seff*altitude);
            case h(3) <= altitude && altitude < h(3+1)
                Seff = S(3)+((altitude-h(3))/100000)*(S(3+1)-S(3));
                rho = rho0*exp(Seff*altitude);
            case h(4) <= altitude && altitude < h(4+1)
                Seff = S(4)+((altitude-h(4))/100000)*(S(4+1)-S(4));
                rho = rho0*exp(Seff*altitude);
            case h(5) <= altitude && altitude < h(5+1)
                Seff = S(5)+((altitude-h(5))/100000)*(S(5+1)-S(5));
                rho = rho0*exp(Seff*altitude); 
            case h(6) <= altitude && altitude < h(6+1)
                Seff = S(6)+((altitude-h(6))/100000)*(S(6+1)-S(6));
                rho = rho0*exp(Seff*altitude);
            case h(7) <= altitude && altitude < h(7+1)
                Seff = S(7)+((altitude-h(7))/100000)*(S(7+1)-S(7));
                rho = rho0*exp(Seff*altitude);
            case h(8) <= altitude && altitude < h(8+1)
                Seff = S(8)+((altitude-h(8))/100000)*(S(8+1)-S(8));
                rho = rho0*exp(Seff*altitude);
            case h(9) <= altitude && altitude < h(9+1)
                Seff = S(9)+((altitude-h(9))/100000)*(S(9+1)-S(9));
                rho = rho0*exp(Seff*altitude);
            case h(10) <= altitude && altitude < h(10+1)
                Seff = S(10)+((altitude-h(10))/100000)*(S(10+1)-S(10));
                rho = rho0*exp(Seff*altitude);
            case altitude >= h(11)
                rho = rho0*exp(S(11)*altitude);                
        end

    elseif contains(atmosphere, "nrlmsise00") 
        %%              
        x3D = r*cos(Phi)*cos(Theta);      % Spherical Coordinates X Coordinate in [m]
        y3D = r*sin(Phi)*cos(Theta);      % Spherical Coordinates Y Coordinate in [m]
        z3D = r*sin(Theta);               % Spherical Coordinates Z Coordinate in [m]  
        P = [x3D y3D z3D];
        
        % Extract time Values
        dt = AtmospericData.DateTime;
        dt.Second = dt.Second + t;
        actdayOfYear = day(dt,'dayofyear');
        actUTseconds = second(dt,'secondofday');
        len = length(dt.Year);
        
        % Calculate DCM ECI to ECEF         
        dcm = dcmeci2ecef('IAU-2000/2006',[dt.Year dt.Month dt.Day dt.Hour dt.Minute dt.Second]);
        % Calculate position in ECEF coordinates
        tmp = arrayfun(@(k) (dcm(:,:,k)*P(k,:)'),1:len,'UniformOutput',false);
        ecef = cell2mat(tmp)';   
        % Calculate LLA coordinates
        lla = ecef2lla(ecef);
        altitude = lla(3);
        
        % Calculate density
        [~, rho] = atmosnrlmsise00(altitude,lla(1),lla(2),dt.Year,actdayOfYear,actUTseconds, AtmospericData.f107a, AtmospericData.f107, ...  
        AtmospericData.aphV);
        rho=rho(6);
                
    end

    % Write solution vector
    rhoV(end+1,1) = rho;
    tV(end+1,1) = t;
    hV(end+1,1) = altitude;

    elseif nargin == 1 % Code for getting solution vector out of ODE   
        rho = rhoV;
        varargout{1} = tV;
        varargout{2} = hV;    
    end
end