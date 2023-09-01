function dfdx = EoM3D (t,x,Ra,Rb,mu,J2,OS,AtmospericData,GravityModel)

    if contains(AtmospericData.atmosphere, "None")
        aD = zeros(1,3);
        aL = zeros(1,3);
    elseif contains(AtmospericData.atmosphere, "Exponential") || contains(AtmospericData.atmosphere, "nrlmsise00")        
        R = sqrt((Ra*cos(x(3)))^2+(Rb*sin(x(3)))^2);            % WGS84 ellipsoid radius [m]
        if isnan(x(1))|| isnan(R)
            StopIn = 1;
            stopode(StopIn)
            rho = -1;
        else   
            rho = Atmosphere(AtmospericData,x(1),x(2),x(3),t,R);    % rho in [kg/m^2]   
        end

        Ds = 86164.0989; % Seconds of siderial day

        aD(1) = -0.5*OS.b^-1*rho*x(4)*sqrt((x(4))^2+(x(1)*cos(x(3))*(x(5)-(2*pi/Ds)))^2+(x(1)*x(6))^2);             % Dragcomponent in r direction[m/s^2]
        aD(2) = -0.5*OS.b^-1*rho*(x(5)-(2*pi/Ds))*sqrt((x(4))^2+(x(1)*cos(x(3))*(x(5)-(2*pi/Ds)))^2+(x(1)*x(6))^2); % Dragcomponent in phi direction [rad/s^2]
        aD(3) = -0.5*OS.b^-1*rho*x(6)*sqrt((x(4))^2+(x(1)*cos(x(3))*(x(5)-(2*pi/Ds)))^2+(x(1)*x(6))^2);             % Dragcomponent in theta direction [rad/s^2]   
        
        
        aL = zeros(1,3);
         

        if OS.LtD ~= 0
            %Method 1
            uL = aD .* OS.LtD;                                      % Adapt Lift with LtD factor
            aL(1) = sqrt((uL(2)*cos(x(3))*x(1))^2+(uL(3)*x(1))^2); % Get radial component through length of theta and phi resulting force
            aL(2) = -(x(1)*cos(x(3))*uL(2)/aL(1))*uL(1);             % Normalize phi and give it length of former radial component 
            aL(3) = -(x(1)*uL(3)/aL(1))*uL(1);                       % Normalize theta and give it length of former radial component 
         
            %Method 2
%             xr(1) = cos(x(2))*cos(x(3));  % Spherical coordinates in [m]
%             xr(2) = sin(x(2))*cos(x(3));  % Spherical coordinates in [m]
%             xr(3) = sin(x(3));            % Spherical coordinates in [m] 
% 
%             v(1) = x(4)*cos(x(2))*cos(x(3))+x(1)*x(5)*cos(x(3))*-sin(x(2))+x(1)*x(6)*-cos(x(2))*sin(x(3));  % Spherical coordinates in [m/s]
%             v(2) = x(4)*sin(x(2))*cos(x(3))+x(1)*x(5)*cos(x(3))*cos(x(2))+x(1)*x(6)*-sin(x(2))*sin(x(3));   % Spherical coordinates in [m/s]
%             v(3) = x(4)*sin(x(3))+x(1)*x(5)*cos(x(3))*0+x(1)*x(6)*cos(x(3));                                % Spherical coordinates in [m/s]
%             vges = sqrt(v(1)^2+v(2)^2+v(3)^2); 
% 
%             vn = v/vges; % normalize velocity vector
% 
%             nO = cross(vn,xr); % finding orbit normal vector
%             eL = cross(nO,vn); % unitvector for Lift
% 
%             aL = eL * OS.LtD * sqrt(aD(1)^2+(aD(2)*x(1))^2+(aD(3)*x(1))^2);
% 
%             rdot0 =aL(1)*cos(x(2))*cos(x(3))+aL(2)*sin(x(2))*cos(x(3))+aL(3)*sin(x(3));
%             Phidot0 =(-aL(1)*sin(x(2))+aL(2)*cos(x(2)));
%             Thetadot0 =(-aL(1)*cos(x(2))*sin(x(3))-aL(2)*sin(x(2))*sin(x(3))+aL(3)*cos(x(3)));
% 
%             aL(1) =rdot0 ;
%             aL(2) =Phidot0 ;
%             aL(3) =Thetadot0 ;
        end
    end
    
    if contains(GravityModel, "Point mass")
        g(1) = -(mu/(x(1)^2));
        g(3) =  0;
    elseif contains(GravityModel, "J2 - Ellipsoid")
        g(1) = -((mu/(x(1)^2))-J2*3*((3*sin(x(3))^2)-1)/((x(1)^4)*2));
        g(3) = -J2*3*cos(x(3))*sin(x(3))/x(1)^5;
    end    

    if sin(x(3))>0.9999 || sin(x(3))<-0.9999
        dfdx=zeros(6,1);
        dfdx(1)= x(4);                                                                          %x1 = r
        dfdx(2)= x(5);                                                                          %x2 = phi
        dfdx(3)= x(6);                                                                          %x3 = theta
        dfdx(4)= x(1)*(x(6)^2)+x(1)*(x(5)^2)*cos(x(3))*cos(x(3))+g(1)+aD(1)+aL(1);              %x4 = rdot
        dfdx(5)= -2*x(4)*x(5)/x(1)+aD(2)+aL(2);                                                 %x5 = phidot
        dfdx(6)= -2*x(4)*x(6)/x(1)-(x(5)^2)*sin(x(3))*cos(x(3))+g(3)+aD(3)+aL(3);               %x6 = thetadot
    else
        dfdx=zeros(6,1);
        dfdx(1)= x(4);                                                                          %x1 = r
        dfdx(2)= x(5);                                                                          %x2 = phi
        dfdx(3)= x(6);                                                                          %x3 = theta
        dfdx(4)= x(1)*(x(6)^2)+x(1)*(x(5)^2)*cos(x(3))*cos(x(3))+g(1)+aD(1)+aL(1);              %x4 = rdot
        dfdx(5)= -2*x(4)*x(5)/x(1)+2*x(5)*x(6)*tan(x(3))+aD(2)+aL(2);                           %x5 = phidot
        dfdx(6)= -2*x(4)*x(6)/x(1)-(x(5)^2)*sin(x(3))*cos(x(3))+g(3)+aD(3)+aL(3);               %x6 = thetadot
    end
    
end