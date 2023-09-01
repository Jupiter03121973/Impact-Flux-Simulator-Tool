function StopOut = stopode(StopIn,varargin)
    persistent stopit;
    if nargin == 0
        StopOut = stopit;
    else
        stopit = StopIn;
    end
end