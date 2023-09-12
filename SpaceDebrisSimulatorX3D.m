function [SimOut,rhoV,aver] = SpaceDebrisSimulatorX3D(AtmospericData,OS,SS,OP,GravityModel,varargin)
    
    global abortSimulation; % Add this line
    abortSimulation = false; % Initialize the flag
    
    clc
    clear Atmosphere
    clear stopode
    
    %% Simulation Parameters
    
    td = SS.td;
    stepsize = SS.stepsize;
    hfinal = SS.hfinal;
    
    %% General Constants
    
    J2 = 1.75553*1e25;  % J2 term for Earth [m^5/s^2]
    mu = 398600.4418*1e9;  % gravitational parameter Earth in [m^3/s^2]
    RE = 6378.14*1e3;    % Earth radius in [m]
    Ra = 6378137.0;     % Earth Semi-major axis in [m] (WGS-84 Ellipsoid)
    Rb = 6356752.3142;  % Earth Semi-minor axis in [m] (WGS-84 Ellipsoid)

    %% Convert Units

    tmax = td*86400;        % max integration time in [s] 
    if stepsize > 0 
        tint = 0:stepsize:tmax; % integration time span and stepsize in [s]
    else
        tint = [0 tmax];        % integration time span (Ode choses own steps) [s]
    end
    
    %% Initial Values
    if OP.RS == 1
                      
        Phi0 = OP.rV(2);                                % [rad]
        Theta0 = OP.rV(3);                              % [rad]
        rdot0 = OP.vV(1);                               % [m]
        Phidot0 = OP.vV(2)*1e-6;                        % [murad/s] to [rad/s]
        Thetadot0 = OP.vV(3)*1e-6;                      % [murad/s] to [rad/s]
        R = sqrt((Ra*cos(Theta0))^2+(Rb*sin(Theta0))^2);% WGS84 Ellipsoid radius [m]
        h0 = (OP.rV(1)*1e3)-R;                          % height [m]
        r0 = R+h0;                                      % initial orbit radius in [m]  
        
    elseif OP.RS == 2
        
        % Orbit Elements to Spherical Coordinates
        OP.a = OP.a*1000;
        
        if OP.e  == 0 || OP.i == 0 % Orbitgeschwindigkeit ist bei berücksichtigung von J2 nicht korrekt
            %[r_ijk,v_ijk] = keplerian2ijk(OP.a,OP.e,OP.i,OP.LAN,OP.AP,OP.f,OP.addN,OP.addV); %alte Funktion
            [r_ijk,v_ijk] = kepler_to_cartesian2(OP.a,OP.e,OP.i,OP.LAN,OP.AP,OP.f, mu, OP.addN,OP.addV);
        else
            %[r_ijk,v_ijk] = keplerian2ijk(OP.a,OP.e,OP.i,OP.LAN,OP.AP,OP.f); %alte Funktion
            [r_ijk,v_ijk] = kepler_to_cartesian(OP.a,OP.e,OP.i,OP.LAN,OP.AP,OP.f, mu);
        end

        % IJK to Spherical Coordinates
        r0 = sqrt(r_ijk(1)^2+r_ijk(2)^2+r_ijk(3)^2);
        Phi0 = atan2(r_ijk(2),r_ijk(1));
        Theta0 = asin(r_ijk(3)/r0);
        rdot0 =v_ijk(1)*cos(Phi0)*cos(Theta0)+v_ijk(2)*sin(Phi0)*cos(Theta0)+v_ijk(3)*sin(Theta0);
        Phidot0 =(-v_ijk(1)*sin(Phi0)+v_ijk(2)*cos(Phi0))/r0;
        Thetadot0 = (v_ijk(3) - rdot0 * cos(Theta0)) / (r0 * cos(Theta0)); % Formel für Thetadot0 korrigiert
        %Thetadot0 =(-v_ijk(1)*cos(Phi0)*sin(Theta0)-v_ijk(2)*sin(Phi0)*sin(Theta0)+v_ijk(3)*cos(Theta0))/r0;     
    end
    
    %% Solve ODEs

    Opt = odeset('Reltol',SS.RelTol,'AbsTol',SS.AbsTol,'Events',@myEvent2,'OutputFcn',@odeprog);
        
    switch SS.ODE
        case '45'
            [t,sc] = ode45(@(t,x) EoM3D(t,x,Ra,Rb,mu,J2,OS,AtmospericData,GravityModel),tint,[r0 Phi0 Theta0 rdot0 Phidot0 Thetadot0],Opt);
        case '23'
            [t,sc] = ode23(@(t,x) EoM3D(t,x,Ra,Rb,mu,J2,OS,AtmospericData,GravityModel),tint,[r0 Phi0 Theta0 rdot0 Phidot0 Thetadot0],Opt);
        case '87'
            [t,sc] = ode87(@(t,x) EoM3D(t,x,Ra,Rb,mu,J2,OS,AtmospericData,GravityModel),tint,[r0 Phi0 Theta0 rdot0 Phidot0 Thetadot0],Opt);
        case '113'
            [t,sc] = ode113(@(t,x) EoM3D(t,x,Ra,Rb,mu,J2,OS,AtmospericData,GravityModel),tint,[r0 Phi0 Theta0 rdot0 Phidot0 Thetadot0],Opt);
        case '15s'
            [t,sc] = ode15s(@(t,x) EoM3D(t,x,Ra,Rb,mu,J2,OS,AtmospericData,GravityModel),tint,[r0 Phi0 Theta0 rdot0 Phidot0 Thetadot0],Opt);
        case '23s'
            [t,sc] = ode23s(@(t,x) EoM3D(t,x,Ra,Rb,mu,J2,OS,AtmospericData,GravityModel),tint,[r0 Phi0 Theta0 rdot0 Phidot0 Thetadot0],Opt);
        case '23t'
            [t,sc] = ode23t(@(t,x) EoM3D(t,x,Ra,Rb,mu,J2,OS,AtmospericData,GravityModel),tint,[r0 Phi0 Theta0 rdot0 Phidot0 Thetadot0],Opt);
        case '23tb'
            [t,sc] = ode23tb(@(t,x) EoM3D(t,x,Ra,Rb,mu,J2,OS,AtmospericData,GravityModel),tint,[r0 Phi0 Theta0 rdot0 Phidot0 Thetadot0],Opt);  
        case '78'
            [t,sc] = ode78(@(t,x) EoM3D(t,x,Ra,Rb,mu,J2,OS,AtmospericData,GravityModel),tint,[r0 Phi0 Theta0 rdot0 Phidot0 Thetadot0],Opt); 
        case '89'
            [t,sc] = ode89(@(t,x) EoM3D(t,x,Ra,Rb,mu,J2,OS,AtmospericData,GravityModel),tint,[r0 Phi0 Theta0 rdot0 Phidot0 Thetadot0],Opt); 
    end
    
    function [value, isterminal, direction] = myEvent2(~,x)
        StopOut = stopode();
        if StopOut == 1
            value = [1;1];
        elseif abortSimulation == 1
            value = [1;1];
        else
            R = sqrt((Ra*cos(x(3)))^2+(Rb*sin(x(3)))^2);
            value=[double(ishandle(95));(x(1)-R <= hfinal*1000)];
        end
        isterminal=[1;1];
        direction=[0;0];
    end

    %% Calculate Additional Values
    
    % Extract Values
    r = sc(:,1);        % [m]
    Phi = sc(:,2);      % [rad]
    Theta = sc(:,3);    % [rad]
    vr = sc(:,4);       % [m/s]
    vphi = sc(:,5);     % [rad/s]
    vtheta = sc(:,6);   % [rad/s]
    
    VgesSp = sqrt(vr.^2+(r.*vphi).^2+(r.*vtheta).^2)/1000; %!!!!!!!!!!!!!!!!!!!!!!
    
    % Get Atmosphere Data
    [rho,tV,hV]= Atmosphere(AtmospericData);
    rhoV = [rho,tV,hV];
    clear Atmosphere  
    
    % 3D Spherical Coordinates 
    x3D = r.*cos(Phi).*cos(Theta);  % Spherical Coordinates Coordinates in [m]
    y3D = r.*sin(Phi).*cos(Theta);  % Spherical Coordinates Coordinates in [m]
    z3D = r.*sin(Theta);            % Spherical Coordinates Coordinates in [m] 
    Vx = vr.*cos(Phi).*cos(Theta)+r.*vphi.*cos(Theta).*-sin(Phi)+r.*vtheta.*-cos(Phi).*sin(Theta);  % Spherical Coordinates Coordinates in [m/s]
    Vy = vr.*sin(Phi).*cos(Theta)+r.*vphi.*cos(Theta).*cos(Phi)+r.*vtheta.*-sin(Phi).*sin(Theta);   % Spherical Coordinates Coordinates in [m/s]
    Vz = vr.*sin(Theta)+r.*vphi.*cos(Theta).*0+r.*vtheta.*cos(Theta);                               % Spherical Coordinates Coordinates in [m/s]
    Vges =sqrt(Vx.^2+Vy.^2+Vz.^2); 
  
    % Orbital specific energie
    if contains(GravityModel, "Point mass")
        epsilon = (0.5.*Vges.^2)-(mu.*r.^-1);                                         % Orbital specific energie
    elseif contains(GravityModel, "J2 - Ellipsoid")
        epsilon = (0.5.*Vges.^2)+J2.*(3.*sin(Theta).^2-1)./(2.*r.^3)-(mu.*r.^-1); % Orbital specific energie with J2
    end
    
    % Orbital parameters
    for in=1:size(x3D)     
        r_ijk = [x3D(in); y3D(in); z3D(in)]; % Spherical Coordinates position vector
        v_ijk = [Vx(in); Vy(in); Vz(in)];    % Spherical Coordinates velocity vector
        
        try
            [~,e,i,LAN,AP,f,TrueLon,arglat,lonper] = ijk2keplerian(r_ijk, v_ijk);    % alternativ function rv2orbelem
            SimOut.a(in)= -mu/(2*epsilon(in));
            %SimOut.a(in) = a;           % Semi-major axis
            SimOut.e(in) = e;           % Orbit eccentricity
            SimOut.i(in) = i;           % Inclination [°]
            SimOut.LAN(in) = LAN;       % Longitude of the ascending node [°]
            SimOut.AP(in) = AP;         % argument of pericenter [°]
            SimOut.f(in) = f;           % Angle between periapsis and current position [°] 
            SimOut.TL(in) = TrueLon;    % Angle between x-axis and position vector [°]
            SimOut.AL(in) = arglat;     % Angle between ascending node and position vector [°]
            SimOut.LP(in) = lonper;     % Angle between x-axis and eccentricity vector [°]
        catch
            SimOut.a(in) = NaN;         % Semi-major axis
            SimOut.e(in) = NaN;         % Orbit eccentricity
            SimOut.i(in) = NaN;         % Inclination [°]
            SimOut.LAN(in) = NaN;       % Longitude of the ascending node [°]
            SimOut.AP(in) = NaN;        % argument of pericenter [°]
            SimOut.f(in) = NaN;         % Angle between periapsis and current position [°] 
            SimOut.TL(in) = NaN;        % Angle between x-axis and position vector [°]
            SimOut.AL(in) = NaN;        % Angle between ascending node and position vector [°]
            SimOut.LP(in) = NaN;        % Angle between x-axis and eccentricity vector [°]
        end
    end    
   
    % 2D Spherical Coordinates
    altitude = r-RE;      % [m] Höhe auf Variablen Wert beziehen!

    if OP.RS == 2
        if OP.i == 0 && OP.e > 0
            x2D = r.*cos((SimOut.LP+SimOut.f).'.*(pi/180)); % [m]
            y2D = r.*sin((SimOut.LP+SimOut.f).'.*(pi/180)); % [m] 
        elseif OP.i == 0 && OP.e == 0
            x2D = r.*cos((SimOut.TL).'.*(pi/180)); % [m]
            y2D = r.*sin((SimOut.TL).'.*(pi/180)); % [m]    
        else
            x2D = r.*cos((SimOut.LAN+SimOut.AP+SimOut.f).'.*(pi/180)); % [m]
            y2D = r.*sin((SimOut.LAN+SimOut.AP+SimOut.f).'.*(pi/180)); % [m]       
        end
    else
        if SimOut.i(1) == 0 && SimOut.e(1) > 0
            x2D = r.*cos((SimOut.LP+SimOut.f).'.*(pi/180)); % [m]
            y2D = r.*sin((SimOut.LP+SimOut.f).'.*(pi/180)); % [m] 
        elseif SimOut.i(1) == 0 && SimOut.e(1) == 0
            x2D = r.*cos((SimOut.TL).'.*(pi/180)); % [m]
            y2D = r.*sin((SimOut.TL).'.*(pi/180)); % [m]    
        else
            x2D = r.*cos((SimOut.LAN+SimOut.AP+SimOut.f).'.*(pi/180)); % [m]
            y2D = r.*sin((SimOut.LAN+SimOut.AP+SimOut.f).'.*(pi/180)); % [m]       
        end
    end

    %%
    % Convert Values   
    r = r/1000;                 % [m] to [km]
    Phi = Phi.*(180/pi);        % [°] to [rad]
    Theta = Theta.*(180/pi);    % [°] to [rad]
    
    SimOut.a= SimOut.a/1000;      % [m] to [km]


    x2D = x2D/1000;             % [m] to [km]
    y2D = y2D/1000;             % [m] to [km]
    
    x3D = x3D/1000;             % [m] to [km]
    y3D = y3D/1000;             % [m] to [km]
    z3D = z3D/1000;             % [m] to [km]
    
    Vges = Vges/1000;           % [m/s] to [km/s]
    t = t/86400;                % [sec] to [d]
    altitude = altitude/1000;   % [m] to [km]
    
    %% filter waves
    
    rows = zeros(1,1);
    trows = zeros(1,1);
    aver = zeros(1,2);
    taver= zeros(1,1);
    n = 1;
    m = 1;
    for in=1:(size(SimOut.f,2)-1)
        if SimOut.f(in) > SimOut.f(in+1)
            rows(n,1) = SimOut.e(in);
            rows(n,2) = r(in);
            trows(n) = t(in);          
            aver(m,:) = mean(rows);
            taver(m,1) = trows(1)+(trows(end)-trows(1))/2;
            rows = zeros(1,1);
            trows = zeros(1,1);
            m = m+1;
            n = 1;
        else
            rows(n,1) = SimOut.e(in); 
            rows(n,2) = r(in);
            trows(n) = t(in);
            n = n+1;
        end
    end
    rows(n,1) = SimOut.e(end);
    rows(n,2) = r(end);
    trows(n) = t(in);
    aver(m,:) = mean(rows);
    taver(m,1) = trows(1)+(trows(end)-trows(1))/2;
    aver = [taver,aver];
    if size(aver,1) > 3
    aver(1,:) = [];
    aver(end,:) = [];
    end
    
    %%
    SimOut.r = r;                % r in [m]
    SimOut.Phi = Phi;            % Phi in [rad]
    SimOut.Theta = Theta;        % Theta in [rad]
    SimOut.x2D = x2D;            % x2D in [m]
    SimOut.y2D = y2D;            % y2D in [m]
    SimOut.x3D = x3D;            % x3D in [km]
    SimOut.y3D = y3D;            % y3D in [km]
    SimOut.z3D = z3D;            % z3D in [km]    
    SimOut.Vges = Vges;          % v ges in [km/s]
    SimOut.VgesSp = VgesSp;
    
    SimOut.t = t;                % t in [d]
    SimOut.altitude = altitude;  % altitude in [km]
    SimOut.epsilon = epsilon;    % specific energie in [m^2/s^2]

end


function [R, V] = kepler_to_cartesian(a, e, i_deg, LAN_deg, AP_deg, f_deg, mu)
    % Convert angles from degrees to radians
    i = deg2rad(i_deg);
    LAN = deg2rad(LAN_deg);
    AP = deg2rad(AP_deg);
    f = deg2rad(f_deg);

    % Improved calculation of specific angular momentum
    %h = sqrt(mu * a * (1 - e) * (1 + e));
    
    % Improved calculation of radius
    r = a * (1 - e^2) / (1 + e * cos(f));
    
    % Velocity vector in the perifocal system
    v_r = sqrt(mu / (a * (1 - e^2))) * e * sin(f);
    v_theta = sqrt(mu / (a * (1 - e^2))) * (1 + e * cos(f));
    
    % Position and velocity vectors in the perifocal system
    r_vec = [r * cos(f); r * sin(f); 0];
    v_vec = [v_r * cos(f) - v_theta * sin(f); v_r * sin(f) + v_theta * cos(f); 0];
    
    % Transformation matrix
    Q = [cos(LAN) * cos(AP) - sin(LAN) * sin(AP) * cos(i), ...
        -cos(LAN) * sin(AP) - sin(LAN) * cos(AP) * cos(i), ...
        sin(LAN) * sin(i);
        
        sin(LAN) * cos(AP) + cos(LAN) * sin(AP) * cos(i), ...
        -sin(LAN) * sin(AP) + cos(LAN) * cos(AP) * cos(i), ...
        -cos(LAN) * sin(i);
        
        sin(AP) * sin(i), cos(AP) * sin(i), cos(i)];
    
    % Transformation of vectors to the inertial frame
    R = Q * r_vec;
    V = Q * v_vec;
end

function [R, V] = kepler_to_cartesian2(a, e, i_deg, LAN_deg, AP_deg, f_deg, mu, addN, addV)
    % Convert angles from degrees to radians
    i = deg2rad(i_deg);
    LAN = deg2rad(LAN_deg);
    AP = deg2rad(AP_deg);
    f = deg2rad(f_deg);

    % Special cases for e and i
    if e == 0 && i ~= 0 && strcmp(addN, 'arglat')
        AP = 0;
        f = addV;  % Argument of Latitude = Argument of Periapsis + True Anomaly
    elseif i == 0 && e ~= 0 && strcmp(addN, 'lonper')
        LAN = 0;
        AP = addV;  % Longitude of Pericenter = Longitude of Ascending Node + Argument of Periapsis
    elseif e == 0 && i == 0 && strcmp(addN, 'TrueLon')
        LAN = 0;
        AP = 0;
        f = addV;  % True Longitude = Longitude of Ascending Node + True Anomaly
    end

    % Improved calculation of specific angular momentum
    %h = sqrt(mu * a * (1 - e) * (1 + e));
    
    % Improved calculation of radius
    r = a * (1 - e^2) / (1 + e * cos(f));
    
    % Velocity vector in the perifocal system
    v_r = sqrt(mu / (a * (1 - e^2))) * e * sin(f);
    v_theta = sqrt(mu / (a * (1 - e^2))) * (1 + e * cos(f));
    
    % Position and velocity vectors in the perifocal system
    r_vec = [r * cos(f); r * sin(f); 0];
    v_vec = [v_r * cos(f) - v_theta * sin(f); v_r * sin(f) + v_theta * cos(f); 0];
    
    % Transformation matrix
    Q = [cos(LAN) * cos(AP) - sin(LAN) * sin(AP) * cos(i), ...
        -cos(LAN) * sin(AP) - sin(LAN) * cos(AP) * cos(i), ...
        sin(LAN) * sin(i);
        
        sin(LAN) * cos(AP) + cos(LAN) * sin(AP) * cos(i), ...
        -sin(LAN) * sin(AP) + cos(LAN) * cos(AP) * cos(i), ...
        -cos(LAN) * sin(i);
        
        sin(AP) * sin(i), cos(AP) * sin(i), cos(i)];
    
    % Transformation of vectors to the inertial frame
    R = Q * r_vec;
    V = Q * v_vec;
end

