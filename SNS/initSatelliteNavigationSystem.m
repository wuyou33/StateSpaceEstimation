function [varargout] = initSatelliteNavigationSystem(func, varargin)
    
    switch func
        case 'init'
            model = init(varargin{1});
            varargout{1} = model;
        otherwise
            error(['Function ''' func ''' not supported.']);
    end
end

function res = init(initArgs)
    res = SatelliteNavigationSystemMock([initArgs.Trajectory initArgs.Velocity]);
end